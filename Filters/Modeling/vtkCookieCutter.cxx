/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCookieCutter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCookieCutter.h"

#include "vtkCellArray.h"
#include "vtkExecutive.h"
#include "vtkGarbageCollector.h"
#include "vtkLine.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkStreamingDemandDrivenPipeline.h"

#include <vector>

vtkStandardNewMacro(vtkCookieCutter);

//----------------------------------------------------------------------------
namespace {

  // Note on the definition of parametric coordinates: Given a sequence of
  // lines segments (vi,vi+1) that form a primitive (e.g., polyline or
  // polygon), the parametric coordinate t along the primitive is
  // [i,i+1). Any point (like an intersection point on the segment) is i+t,
  // where 0 <= t < 1.


  // Infrastructure for cropping----------------------------------------------
  struct SortPoint
  {
    double Tp; //parametric coordinate along primitive (line or polygon)
    double Tl; //parametric coordinate along loop (the cookie cutter)
    double X[3]; //position of point
    int Type; // type of vertex (intersection, vertex)
    int Class; // inside/outside classification of vertex
    enum TypeEnum {INTERSECTION=0,VERTEX=1}; //type of point
    enum ClassEnum {UNKNOWN=0, INSIDE=1,OUTSIDE=2}; //classification of segment (i,(i+1)%npts)
    SortPoint(double tp, double tl, double x[3], int type) :
      Tp(tp), Tl(tl), Type(type)
    {
      X[0]=x[0]; X[1]=x[1]; X[2]=x[2];
      this->Class = SortPoint::UNKNOWN;
    }
  };

  // Special sort operation on primitive parametric coordinate Tp------------
  bool PointSorter(SortPoint const& lhs, SortPoint const& rhs)
  {
    return lhs.Tp < rhs.Tp;
  }

  // Vectors are used to hold points.
  typedef std::vector<SortPoint> SortPointType;

  // Clean up list of intersections -----------------------------------------
  // The points are assumed sorted in parametric coordinates, not closed
  int CleanSortedPolyline(SortPointType &sortedPoints)
  {
    int i, j;
    double ti, tj, delT;
    for ( i=0; i < sortedPoints.size()-1; )
    {
      j = i + 1;
      ti = sortedPoints[i].Tp;
      tj = sortedPoints[j].Tp;

      delT = tj - static_cast<double>(static_cast<int>(ti));
      if ( delT <= 0.001 ) //arbitrary epsilon
      {
        sortedPoints.erase(sortedPoints.begin()+j);
        sortedPoints[i].Type = SortPoint::INTERSECTION;
      }
      else
      {
        ++i;
      }
    }
    return static_cast<int>(sortedPoints.size());
  }

  // Clean up list of intersections -----------------------------------------
  // The points are assumed sorted in parametric coordinates, closed loop
  int CleanSortedPolygon(vtkIdType npts, SortPointType &sortedPoints)
  {
    if ( sortedPoints.size() < 2 )
    {
      return static_cast<int>(sortedPoints.size());
    }

    int i, j;
    double ti, tj;
    size_t num;
    for ( i=0; i < sortedPoints.size(); )
    {
      num = sortedPoints.size();
      j = (i + 1) % num;
      ti = sortedPoints[i].Tp;
      tj = sortedPoints[j].Tp;
      tj = (tj < ti ? tj += static_cast<double>(npts) : tj);

      if ( (tj-ti) <= 0.001 ) //arbitrary epsilon
      {
        sortedPoints.erase(sortedPoints.begin()+j);
      }
      else
      {
        ++i;
      }
    }
    return static_cast<int>(sortedPoints.size());
  }

  // Convenience method------------------------------------------------------
  inline void InsertPoint(double x[3], vtkPoints *outPts, vtkCellArray *ca,
                          vtkCellData *inCD, vtkCellData *outCD,
                          vtkIdType cellId, vtkIdType newCellId)
  {
    vtkIdType newPtId = outPts->InsertNextPoint(x);
    ca->InsertCellPoint(newPtId);
    outCD->CopyData(inCD,cellId,newCellId);
  }

  // What type of segment are we processing
  enum SegmentType {Poly, Loop};

  // Compute interval-------------------------------------------------------
  // Compute the number of points between two parametric coordinates in
  // modulo coordinates around a polygon, and provide the information
  // necessary to extract points. The trick is going in the right
  // direction. The "orient" variableindicates whether the loop and polygon
  // are oriented in the same way which is important for certina operations.
  int ComputeTRange(double orient, int npts, int num, SortPointType &sortedPoints,
                    vtkIdType spi, vtkIdType spj, SegmentType type, int &idx, int &inc)
  {
    int i, j;

    // We alsways go around the polygon in the positive direction. Hence we
    // look at the parametric coordinates Tp.
    if ( type == SegmentType::Poly )
    {
      i = static_cast<int>(sortedPoints[spi].Tp);
      j = static_cast<int>(sortedPoints[spj].Tp);

      // always move in the forward directions
      inc = 1;
      idx = (i+1) % npts;
      j = ( j < i ? j+npts : j ); //modulo forward
      if ( num != 2 || i != j )
      {
        return (j - i);
      }
      // A specical case exists when it's not clear whether the return around
      // the loop is along the (i,j) edge or around the cutting loop. This
      // only occurs when the loop cuts the same polygon segment(i==j) in
      // two places.
      else
      {
        if ( orient > 0.0 )
        {
          if ( sortedPoints[0].Tl < sortedPoints[1].Tl )
          {
            idx = (i + 1) % npts;
            return npts;
          }
          else
          {
            return 0;
          }
        }
        else //loops are in opposite directions
        {
          if ( sortedPoints[0].Tl > sortedPoints[1].Tl )
          {
            idx = (i + 1) % npts;
            return npts;
          }
          else
          {
            return 0;
          }
        }
      }
    }//if processing the polygon

    // Else if processing the cutting loop, use Tl
    else // type == SegmentType::Loop
    {
      i = static_cast<int>(sortedPoints[spi].Tl);
      j = static_cast<int>(sortedPoints[spj].Tl);
      if ( orient > 0.0 ) //forward direction
      {
        inc = 1;
        idx = (i+1) % npts;
        j = ( j < i ? j+npts : j ); //modulo forward
        return (j - i);
      }
      else //reverse direction
      {
        inc = (-1);
        idx = i;
        i = ( i < j ? i+npts : i );
        return (i - j);
      }
    }//around cutting loop
  }

  // Process a polyline-------------------------------------------------------
  void CropLine(vtkIdType cellId, vtkIdType npts, vtkIdType *pts, vtkPoints *inPts,
                vtkPolygon *loop, double *p, double bds[6], double n[3],
                vtkCellData *inCellData, vtkPoints *outPts,
                vtkCellArray *outLines, vtkCellData *outCellData)
  {
    // Make sure that this is a valid line
    if ( npts < 2 )
    {
      return;
    }

    // Create a vector of points with parametric coordinates, etc. which will
    // be sorted later. The tricky part is getting the classification of each
    // point correctly.
    // First insert all of the points defining the polyline
    vtkIdType i, j, newCellId;
    double t, u, v, x[3], x0[3], x1[3], y0[3], y1[3];
    SortPointType sortedPoints;
    for (i=0; i < npts; ++i)
    {
      t = static_cast<double>(i);
      inPts->GetPoint(pts[i],x);
      sortedPoints.push_back(SortPoint(t,0.0,x,SortPoint::VERTEX));
    }

    // Now insert any intersection points
    vtkIdType numInts=0, numLoopPts=loop->Points->GetNumberOfPoints();
    for (numInts=0, i=0; i < (npts-1); ++i)
    {
      inPts->GetPoint(pts[i],x0);
      inPts->GetPoint(pts[i+1],x1);

      // Traverse polygon loop intersecting each polygon segment
      for (j=0; j < numLoopPts; ++j)
      {
        loop->Points->GetPoint(j,y0);
        loop->Points->GetPoint((j+1)%numLoopPts,y1);
        if ( vtkLine::Intersection(x0,x1,y0,y1,u,v) == 2 )
        {
          numInts++;
          x[0] = x0[0] + u*(x1[0]-x0[0]);
          x[1] = x0[1] + u*(x1[1]-x0[1]);
          x[2] = x0[2] + u*(x1[2]-x0[2]);
          u += static_cast<double>(i);
          v += static_cast<double>(j);
          sortedPoints.push_back(SortPoint(u,v,x,SortPoint::INTERSECTION));
        }
      }//intersect all line segments that form the loop
    }//for all line segments that make up this polyline

    // Okay sort that sucker in parametric coordinates
    vtkIdType num = sortedPoints.size();
    std::sort(sortedPoints.begin(), sortedPoints.end(), &PointSorter);

    // Clean up sorted array, remove nearly conincident points
    if ( numInts > 0 )
    {
      num = CleanSortedPolyline(sortedPoints);
    }

    // Determine the classification of the first segment of the
    // polyline. This will set us up for later processing.
    double midSegment[3];
    midSegment[0] = 0.5 * (sortedPoints[0].X[0] + sortedPoints[1].X[0]);
    midSegment[1] = 0.5 * (sortedPoints[0].X[1] + sortedPoints[1].X[1]);
    midSegment[2] = 0.5 * (sortedPoints[0].X[2] + sortedPoints[1].X[2]);

    if ( vtkPolygon::PointInPolygon(midSegment, numLoopPts, p, bds, n) == 1 )
    {
      sortedPoints[0].Class = SortPoint::INSIDE;
    }
    else
    {
      sortedPoints[0].Class = SortPoint::OUTSIDE;
    }

    // If no intersection points, then the line is either entirely in or
    // entirely out
    if ( numInts == 0 )
    {
      if ( sortedPoints[0].Class == SortPoint::OUTSIDE )
      {
        return;
      }
      else //whole line is inside and therefore output
      {
        newCellId = outLines->InsertNextCell(npts);
        for (i=0; i < npts; ++i)
        {
          InsertPoint(inPts->GetPoint(pts[i]), outPts, outLines,
                      inCellData, outCellData, cellId, newCellId);
        }
      }
    }//if no intersections

    // There is at least one intersection. Propagate the initial
    // classification along the polyline.
    int currentClass = sortedPoints[0].Class;
    for (i=1; i < num; ++i)
    {
      if ( sortedPoints[i].Type == SortPoint::INTERSECTION )
      {
        // Flip classification
        currentClass = (currentClass == SortPoint::INSIDE ?
                        SortPoint::OUTSIDE : SortPoint::INSIDE);
      }
      sortedPoints[i].Class = currentClass;
    }

    // If here, then pieces of the intersected polyline are output.
    vtkIdType startIdx=0, endIdx=0, numInsertedPts;
    while ( startIdx < (num-1) )
    {
      // move to the beginning of the next interior segment.
      while ( startIdx < (num-1) && sortedPoints[startIdx].Class == SortPoint::OUTSIDE )
      {
        startIdx++;
      }
      if ( startIdx >= (num-1) )
      {
        continue; // all done
      }

      // Find the end of the run, i.e., link together interior segments.
      endIdx = startIdx + 1;
      while ( endIdx < (num-1) && sortedPoints[endIdx].Class == SortPoint::INSIDE )
      {
        endIdx++;
      }

      // Output line segments
      numInsertedPts = endIdx - startIdx + 1;
      newCellId = outLines->InsertNextCell(numInsertedPts);
      for (i=startIdx; i <= endIdx; ++i)
      {
        InsertPoint(sortedPoints[i].X, outPts, outLines,
                    inCellData,outCellData,cellId,newCellId);
      }
      startIdx = endIdx;
    }//over all sorted points
  }//CropLine


  // Process a polyon--------------------------------------------------------
  void CropPoly(vtkIdType cellId, vtkIdType npts, vtkIdType *pts, vtkPoints *inPts,
                vtkPolygon *loop, double *p, double bds[6], double n[3],
                vtkCellData *inCellData, vtkPoints *outPts,
                vtkCellArray *outPolys, vtkCellData *outCellData)
  {
    // Make sure that this is a valid polygon
    if ( npts < 3 )
    {
      return;
    }

    // Run around the polygon inserting intersection points.
    vtkIdType i, j, newCellId;
    double u, v, x[3], x0[3], x1[3], y0[3], y1[3];
    SortPointType sortedPoints;
    vtkIdType numLoopPts=loop->Points->GetNumberOfPoints();
    for (i=0; i < npts; ++i)
    {
      inPts->GetPoint(pts[i],x0);
      inPts->GetPoint(pts[(i+1)%npts],x1);

      // Traverse polygon loop intersecting each polygon segment
      for (j=0; j < numLoopPts; ++j)
      {
        loop->Points->GetPoint(j,y0);
        loop->Points->GetPoint((j+1)%numLoopPts,y1);
        if ( vtkLine::Intersection(x0,x1,y0,y1,u,v) == 2 )
        {
          x[0] = x0[0] + u*(x1[0]-x0[0]);
          x[1] = x0[1] + u*(x1[1]-x0[1]);
          x[2] = x0[2] + u*(x1[2]-x0[2]);
          u += static_cast<double>(i);
          v += static_cast<double>(j);
          sortedPoints.push_back(SortPoint(u,v,x,SortPoint::INTERSECTION));
        }
      }//intersect all line segments that form the loop
    }//for all line segments that make up this polygon

    // The number of intersection points. Sort around the polygon.
    vtkIdType num = sortedPoints.size();
    std::sort(sortedPoints.begin(), sortedPoints.end(), &PointSorter);

    // If there are no intersection points, then the polygon is either
    // entirely in or entirely out. We add the polygon (or not) as
    // appropriate using an in/out test on one of the polygon vertices.
    if ( num <= 0 )
    {
      if ( vtkPolygon::PointInPolygon(inPts->GetPoint(pts[0]), numLoopPts, p, bds, n) == 1 )
      {
        newCellId = outPolys->InsertNextCell(npts);
        for (i=0; i < npts; ++i)
        {
          InsertPoint(inPts->GetPoint(pts[i]), outPts, outPolys,
                      inCellData, outCellData, cellId, newCellId);
        }
      }
      return;
    }//if no intersections
    else
    {
      num = CleanSortedPolygon(npts,sortedPoints);
    }

    // If a weird degenerate case occurs where the number of intersections is
    // not even, then the polygon can be nestled into the "corner" of the
    // loop, i.e., have one of more vertices that lies on a loop edge or
    // vertex.  Do a dirty check by performing an inside check of the center
    // point of the polygon.
    if ( (num % 2) ) // generally there are multiples of two intersections
    {
      double center[3]={0,0,0};
      for ( i=0; i < npts; ++i )
      {
        inPts->GetPoint(pts[i],x);
        center[0] += x[0];
        center[1] += x[1];
        center[2] += x[2];
      }
      center[0] /= static_cast<double>(npts);
      center[1] /= static_cast<double>(npts);
      center[2] /= static_cast<double>(npts);

      if ( vtkPolygon::PointInPolygon(center, numLoopPts, p, bds, n) == 1 )
      {
        newCellId = outPolys->InsertNextCell(npts);
        for (i=0; i < npts; ++i)
        {
          InsertPoint(inPts->GetPoint(pts[i]), outPts, outPolys,
                      inCellData, outCellData, cellId, newCellId);
        }
      }
      return;
    }//if degenerate case

    // If here, there are at least two intersection points. We need to add
    // additional "interior" vertices (non intersection point vertices) that
    // come from either the loop or the polygon we are cutting. We use a
    // simple winding rule: at each intersection point there is a choice as
    // to which path to take. We take the path that turns towards the
    // interior of the polygon (as determined from the cross product). This
    // works becase the loop-poly operation is an intersection set
    // operation. It is much faster and more reliable than PointInPolygon()
    // tests.
    //
    // Currently this code extract one loop (polygon) as a result of the
    // intersection. TODO in the future modify the code to extract multiple
    // polygons.
    int uMin, uMax, vMin, vMax, numInsertedPoints=0;
    double orient, pNormal[3], vx[3], vy[3], vxy[3];
    vtkPolygon::ComputeNormal(inPts,npts,pts,pNormal);

    // Orientation of the two loops with respect to one another.
    orient = ( vtkMath::Dot(n,pNormal) > 0.0 ? 1.0 : -1.0 );

    // Begin the process of inserting cell points. Initially insert the
    // number of intersections and then update it later.
    newCellId = outPolys->InsertNextCell(num);
    int idx, inc, numToInsert, ii;
    for ( i=0; i < num; ++i )
    {
      // Start the current loop here; grab the next interection point.
      j = (i+1) % num;

      // Insert this intersection point
      numInsertedPoints++;
      InsertPoint(sortedPoints[i].X, outPts, outPolys,
                  inCellData, outCellData, cellId, newCellId);

      // Determine which path to take. First grab two line segments
      // that intersect.
      uMin = vtkMath::Floor(sortedPoints[i].Tp);
      uMax = vtkMath::Ceil(sortedPoints[i].Tp) % npts;
      inPts->GetPoint(pts[uMin],x0);
      inPts->GetPoint(pts[uMax],x1);
      vx[0] = x1[0] - x0[0];
      vx[1] = x1[1] - x0[1];
      vx[2] = x1[2] - x0[2];

      vMin = vtkMath::Floor(sortedPoints[i].Tl);
      vMax = vtkMath::Ceil(sortedPoints[i].Tl) % numLoopPts;
      loop->Points->GetPoint(vMin,y0);
      loop->Points->GetPoint(vMax,y1);
      vy[0] = y1[0] - y0[0];
      vy[1] = y1[1] - y0[1];
      vy[2] = y1[2] - y0[2];

      // Choose which segment to traverse, and insert points interior
      // to the segment
      vtkMath::Cross(vx,vy,vxy);
      if ( (orient*vtkMath::Dot(vxy,pNormal)) <= 0 )
      {
        numToInsert = ComputeTRange(orient, npts, num, sortedPoints, i, j,
                                    SegmentType::Poly, idx, inc);
        for ( ii=0; ii < numToInsert; ++ii )
        {
          numInsertedPoints++;
          InsertPoint(inPts->GetPoint(pts[idx]), outPts, outPolys,
                      inCellData, outCellData, cellId, newCellId);
          idx = (idx+inc) % npts; //always going in forward direction
        }
      }
      else
      {
        numToInsert = ComputeTRange(orient, numLoopPts, num, sortedPoints,
                                    i, j, SegmentType::Loop, idx, inc);
        for ( ii=0; ii < numToInsert; ++ii )
        {
          numInsertedPoints++;
          InsertPoint(loop->Points->GetPoint(idx), outPts, outPolys,
                      inCellData, outCellData, cellId, newCellId);
          idx = (idx+inc) % numLoopPts;
          idx = (idx < 0 ? (numLoopPts-idx) : idx);//needed if going backwards
        }
      }

    }//over all intersection points
    outPolys->UpdateCellCount(numInsertedPoints);

  }//CropPoly

} //anonymous namespace

//----------------------------------------------------------------------------
// Instantiate object with empty loop.
vtkCookieCutter::vtkCookieCutter()
{
  this->SetNumberOfInputPorts(2);
}

//----------------------------------------------------------------------------
vtkCookieCutter::~vtkCookieCutter()
{
}

//----------------------------------------------------------------------------
void vtkCookieCutter::SetLoopsConnection(vtkAlgorithmOutput* algOutput)
{
  this->SetInputConnection(1, algOutput);
}

//----------------------------------------------------------------------------
vtkAlgorithmOutput *vtkCookieCutter::GetLoopsConnection()
{
  return this->GetInputConnection(1, 0);
}

//----------------------------------------------------------------------------
void vtkCookieCutter::SetLoopsData(vtkDataObject *input)
{
  this->SetInputData(1, input);
}

//----------------------------------------------------------------------------
vtkDataObject *vtkCookieCutter::GetLoops()
{
  if (this->GetNumberOfInputConnections(1) < 1)
  {
    return NULL;
  }

  return this->GetExecutive()->GetInputData(1, 0);
}

//----------------------------------------------------------------------------
int vtkCookieCutter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *loopInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *loops = vtkPolyData::SafeDownCast(
    loopInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPolyData *output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  // Initialize and check data
  vtkDebugMacro(<<"Cookie cutting...");

  vtkIdType numPts, numLoopPts;
  if ( (numPts = input->GetNumberOfPoints()) < 1 )
  {
    vtkErrorMacro("Input contains no points");
    return 1;
  }

  if ( !loops || loops->GetNumberOfCells() < 1 )
  {
    vtkErrorMacro("Please define a polygonal loop with at least three points");
    return 1;
  }

  // Create output data objects and prepare for processing
  vtkPoints *inPts = input->GetPoints();
  vtkPoints *loopPts = loops->GetPoints();
  vtkPoints *outPts = inPts->NewInstance();

  vtkCellData *inCellData = input->GetCellData();
  vtkCellData *outCellData = output->GetCellData();
  outCellData->CopyAllocate(inCellData);

  vtkCellArray *inVerts = input->GetVerts();
  vtkCellArray *outVerts = vtkCellArray::New();
  outVerts->Allocate(inVerts->GetNumberOfCells(),1);

  vtkCellArray *inLines = input->GetLines();
  vtkCellArray *outLines = vtkCellArray::New();
  outLines->Allocate(inLines->GetNumberOfCells(),1);

  vtkCellArray *inPolys = input->GetPolys();
  vtkCellArray *inStrips = input->GetStrips();
  vtkCellArray *outPolys = vtkCellArray::New();
  outPolys->Allocate(inPolys->GetNumberOfCells(),1);

  // Initialize and create polygon representing the loop
  double n[3], bds[6];
  vtkPolygon *poly = vtkPolygon::New();

  // Process loops from second input. Note that the cell id in vtkPolyData
  // starts with verts, then lines, then polys, then strips.
  vtkIdType i, npts, *pts, cellId=0, newPtId, newCellId;
  vtkCellArray *loopPolys = loops->GetPolys();
  for (loopPolys->InitTraversal(); loopPolys->GetNextCell(npts,pts); )
  {
    numLoopPts = npts;
    if ( numLoopPts < 3 )
    {
      continue; // need a valid polygon, skip this one
    }
    poly->Initialize(npts,pts,loopPts);
    poly->GetBounds(bds);
    vtkPolygon::ComputeNormal(poly->Points,n);
    double *p = static_cast<double*>(poly->Points->GetVoidPointer(0));

    // Start by processing the verts. A simple in/out check.
    double x[3];
    if ( inVerts->GetNumberOfCells() > 0 )
    {
      for (inVerts->InitTraversal(); inVerts->GetNextCell(npts,pts); ++cellId)
      {
        for ( i=0; i<npts; ++i)
        {
          inPts->GetPoint(pts[i], x);
          if ( vtkPolygon::PointInPolygon(x, numLoopPts, p, bds, n) == 1 )
          {
            newPtId = outPts->InsertNextPoint(x);
            newCellId = outVerts->InsertNextCell(1,&newPtId);
            outCellData->CopyData(inCellData,cellId,newCellId);
          }
        }
      }//for all verts
    }//if vert cells

    // Now process lines
    if ( inLines->GetNumberOfCells() > 0 )
    {
      for (inLines->InitTraversal(); inLines->GetNextCell(npts,pts); ++cellId)
      {
        CropLine(cellId, npts,pts,inPts,
                 poly,p,bds,n,inCellData,
                 outPts,outLines,outCellData);
      }
    }//if line cells

    // Now process polygons
    if ( inPolys->GetNumberOfCells() > 0 )
    {
      for (inPolys->InitTraversal(); inPolys->GetNextCell(npts,pts); ++cellId)
      {
        CropPoly(cellId, npts,pts,inPts,
                 poly,p,bds,n,inCellData,
                 outPts,outPolys,outCellData);
      }
    }//if polygonal cells

    // Now process triangle strips
    if ( inStrips->GetNumberOfCells() > 0 )
    {
      for (inStrips->InitTraversal(); inStrips->GetNextCell(npts,pts); ++cellId)
      {
        vtkIdType numTriPts=3, triPts[3];
        for ( i=0; i<(npts-2); ++i)
        {
          triPts[0] = pts[i];
          triPts[1] = pts[i+1];
          triPts[2] = pts[i+2];
          CropPoly(cellId, numTriPts,triPts,inPts,
                   poly,p,bds,n,inCellData,
                   outPts,outPolys,outCellData);
        }
      }
    }//if polygonal cells
  }//for all loops

  // Assign output as appropriate
  output->SetVerts(outVerts);
  outVerts->Delete();

  output->SetLines(outLines);
  outLines->Delete();

  output->SetPolys(outPolys);
  outPolys->Delete();

  // Clean up
  output->SetPoints(outPts);
  outPts->Delete();
  poly->Delete();

  return 1;
}

//----------------------------------------------------------------------------
int vtkCookieCutter::
RequestUpdateExtent(vtkInformation *vtkNotUsed(request),
                    vtkInformationVector **inputVector,
                    vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *loopInfo = inputVector[1]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  if (loopInfo)
  {
    loopInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
                    0);
    loopInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
                    1);
    loopInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
                    0);
  }
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER(),
              outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_PIECE_NUMBER()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES(),
              outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_PIECES()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS(),
              outInfo->Get(vtkStreamingDemandDrivenPipeline::UPDATE_NUMBER_OF_GHOST_LEVELS()));
  inInfo->Set(vtkStreamingDemandDrivenPipeline::EXACT_EXTENT(), 1);

  return 1;
}

//----------------------------------------------------------------------------
int vtkCookieCutter::FillInputPortInformation(int port, vtkInformation *info)
{
  if (port == 0)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  else if (port == 1)
  {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkPolyData");
    return 1;
  }
  return 0;
}

//----------------------------------------------------------------------------
void vtkCookieCutter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

}
