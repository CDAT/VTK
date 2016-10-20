/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCookieCutFilter.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#include "vtkCookieCutFilter.h"

#include "vtkCellArray.h"
#include "vtkExecutive.h"
#include "vtkGarbageCollector.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkLine.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkCellData.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolygon.h"

#include <vector>

vtkStandardNewMacro(vtkCookieCutFilter);
vtkCxxSetObjectMacro(vtkCookieCutFilter,Loop,vtkPoints);

//----------------------------------------------------------------------------
namespace {

  // Infrastructure for interaction
  struct SortPoint
  {
    double T; //parametric coordinate
    double X[3]; // physical location
    int Type; //classification (boundary, vertex)
    enum TypeEnum {BOUNDARY=0,VERTEX=1};
    SortPoint(double t, double x[3], int type) :
      T(t), Type(type)
    {
      X[0]=x[0]; X[1]=x[1]; X[2]=x[2];
    }
  };

  // Special sort operation on parametric coordinate
  bool PointSorter(SortPoint const& lhs, SortPoint const& rhs)
  {
    return lhs.T < rhs.T;
  }

  // Vectors are used to hold points.
  typedef std::vector<SortPoint> SortPointType;

  // Process a polyline
  void CropLine(vtkIdType cellId, vtkIdType npts, vtkIdType *pts, vtkPoints *inPts,
                vtkPolygon *poly, double *p, double bds[6], double n[3],
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
    vtkIdType i, j, newPtId, newCellId;
    double t, u, v, x[3], x0[3], x1[3], y0[3], y1[3];
    SortPointType sortedPoints;
    for (i=0; i < npts; ++i)
    {
      t = static_cast<double>(i);
      inPts->GetPoint(pts[i],x);
      sortedPoints.push_back(SortPoint(t,x,SortPoint::VERTEX));
    }

    // Now insert any intersection points
    vtkIdType numInts=0, numLoopPts=poly->Points->GetNumberOfPoints();
    for (numInts=0, i=0; i < (npts-1); ++i)
    {
      inPts->GetPoint(pts[i],x0);
      inPts->GetPoint(pts[i+1],x1);

      // Traverse polygon loop intersecting each polygon segment
      for (j=0; j < numLoopPts; ++j)
      {
        poly->Points->GetPoint(j,y0);
        poly->Points->GetPoint((j+1)%numLoopPts,y1);
        if ( vtkLine::Intersection(x0,x1,y0,y1,u,v) == 2 )
        {
          numInts++;
          t = static_cast<double>(i) + u;
          sortedPoints.push_back(SortPoint(t,x,SortPoint::BOUNDARY));
        }
      }//intersect all line segments that form the loop
    }//for all line segments that make up this polyline

    // Okay sort that sucker and determine the classification of the
    // first point. This will set us up for later processing.
    vtkIdType num = sortedPoints.size();
    std::sort(sortedPoints.begin(), sortedPoints.end(), &PointSorter);
    if ( vtkPolygon::PointInPolygon(sortedPoints[0].X, numLoopPts, p, bds, n) == 1 )
    {
      sortedPoints[0].Type = SortPoint::BOUNDARY; //beginning of inside run
    }
    sortedPoints[num-1].Type = SortPoint::BOUNDARY; //end of last run

    // If no intersection points, then the line is either entirely in or
    // entirely out
    if ( numInts == 0 )
    {
      if ( sortedPoints[0].Type == SortPoint::VERTEX )
      {
        return;
      }
      else //whole line is output
      {
        newCellId = outLines->InsertNextCell(npts);
        for (i=0; i < npts; ++i)
        {
          newPtId = outPts->InsertNextPoint(inPts->GetPoint(pts[i]));
          outLines->InsertCellPoint(newPtId);
          outCellData->CopyData(inCellData,cellId,newCellId);
        }
      }
    }//no intersections

    // Perform a little sanity check, If something is amiss, don't spit out
    // anything.
    int startType = sortedPoints[0].Type;
    int endType = SortPoint::VERTEX;
    if ( vtkPolygon::PointInPolygon(sortedPoints[num-1].X, numLoopPts, p, bds, n) == 1 )
    {
      endType = SortPoint::BOUNDARY; //last point is part of inside run
    }
    if ( (!(numInts % 2) && endType != startType) ||
         ((numInts % 2) && endType == startType) )
    {
      return; //topological inconsistency
    }

    //If here, then pieces of the intersected line are spit out. Also it means at
    //least one intersection.
    vtkIdType startIdx=0, endIdx=0, numInsertedPts;
    while ( startIdx < (num-1) )
    {
      //move to beginning of valid segment. Remember first point has been classified
      //as either on boundary of valid segment, or part of an outside run.
      while ( sortedPoints[startIdx].Type == SortPoint::VERTEX )
      {
        startIdx++;
      }

      // Find the end of the run, look for next boundary
      endIdx = startIdx + 1;
      while ( endIdx < num && sortedPoints[endIdx].Type == SortPoint::VERTEX )
      {
        endIdx++;
      }

      //spit out segment
      numInsertedPts = endIdx - startIdx + 1;
      newCellId = outLines->InsertNextCell(numInsertedPts);
      for (i=startIdx; i <= endIdx; ++i)
      {
        newPtId = outPts->InsertNextPoint(sortedPoints[i].X);
        outLines->InsertCellPoint(newPtId);
        outCellData->CopyData(inCellData,cellId,newCellId);
      }
      startIdx = endIdx;
    }//over all sorted points
  }//CropLine

} //anonymous namespace

//----------------------------------------------------------------------------
// Instantiate object with empty loop.
vtkCookieCutFilter::vtkCookieCutFilter()
{
  this->Loop = NULL;
}

//----------------------------------------------------------------------------
vtkCookieCutFilter::~vtkCookieCutFilter()
{
  if (this->Loop)
  {
    this->Loop->Delete();
    this->Loop = NULL;
  }
}

//----------------------------------------------------------------------------
int vtkCookieCutFilter::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
  // get the info objects
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  // get the input and output
  vtkPolyData *input = vtkPolyData::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));
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

  if ( this->Loop == NULL ||
  (numLoopPts=this->Loop->GetNumberOfPoints()) < 3 )
  {
    vtkErrorMacro("Please define a loop with at least three points");
    return 1;
  }

  // Initialize and create polygon representing the loop
  vtkPolygon *poly = vtkPolygon::New();
  poly->Initialize(numLoopPts, this->Loop);
  double n[3]={0,0,1}, bds[6];
  poly->GetBounds(bds);
  double *p = static_cast<double*>(poly->Points->GetVoidPointer(0));

  vtkPoints *inPts = input->GetPoints();
  vtkPoints *outPts = inPts->NewInstance();

  vtkCellData *inCellData = input->GetCellData();
  vtkCellData *outCellData = output->GetCellData();
  outCellData->CopyAllocate(inCellData);

  // Start by processing the verts. A simple in/out check.
  vtkIdType i, npts, *pts, cellId, newPtId, newCellId;
  double x[3];
  vtkCellArray *inVerts = input->GetVerts();
  if ( inVerts->GetNumberOfCells() > 0 )
  {
    vtkCellArray *outVerts = vtkCellArray::New();
    outVerts->Allocate(inVerts->GetNumberOfCells(),1);
    for (cellId=0, inVerts->InitTraversal();
         inVerts->GetNextCell(npts,pts); ++cellId)
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
    output->SetVerts(outVerts);
    outVerts->Delete();
  }//if vert cells

  // Now process lines
  vtkCellArray *inLines = input->GetLines();
  if ( inLines->GetNumberOfCells() > 0 )
  {
    vtkCellArray *outLines = vtkCellArray::New();
    outLines->Allocate(inLines->GetNumberOfCells(),1);
    for (cellId=0, inLines->InitTraversal();
         inLines->GetNextCell(npts,pts); ++cellId)
    {
      CropLine(cellId, npts,pts,inPts,
               poly,p,bds,n,inCellData,
               outPts,outLines,outCellData);
    }
    output->SetLines(outLines);
    outLines->Delete();
  }//if line cells

  // Now process polygons
  vtkCellArray *inPolys = input->GetPolys();
  if ( inPolys->GetNumberOfCells() > 0 )
  {
    vtkCellArray *outPolys = vtkCellArray::New();
    outPolys->Allocate(inPolys->GetNumberOfCells(),1);
    for (cellId=0, inPolys->InitTraversal();
         inPolys->GetNextCell(npts,pts); ++cellId)
    {
    }
    output->SetPolys(outPolys);
    outPolys->Delete();
  }//if line cells

  // Clean up
  output->SetPoints(outPts);
  outPts->Delete();
  poly->Delete();

  return 1;
}

//----------------------------------------------------------------------------
vtkMTimeType vtkCookieCutFilter::GetMTime()
{
  vtkMTimeType mTime=this->Superclass::GetMTime();
  vtkMTimeType time;

  if ( this->Loop != NULL )
  {
    time = this->Loop->GetMTime();
    mTime = ( time > mTime ? time : mTime );
  }

  return mTime;
}

//----------------------------------------------------------------------------
void vtkCookieCutFilter::PrintSelf(ostream& os, vtkIndent indent)
{
  this->Superclass::PrintSelf(os,indent);

  if ( this->Loop )
  {
    os << indent << "Loop of " << this->Loop->GetNumberOfPoints()
       << "points defined\n";
  }
  else
  {
    os << indent << "Loop not defined\n";
  }
}
