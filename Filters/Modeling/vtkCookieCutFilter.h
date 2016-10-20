/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCookieCutFilter.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
/**
 * @class   vtkCookieCutFilter
 * @brief   cut vtkPolyData defined on 2D plane with a polyon
 *
 * This filter crops an input vtkPolyData consisting of simplicial cells
 * (i.e., points, lines, and triangles) with a loop specified by an ordered
 * list of points defining a single polygon. The output consists of the
 * cropped simplicial cells.
 *
 * @warning
 * The z-values of the input vtkPolyData and the points defining the loop are
 * assumed to lie at z=0. In other words, this filter assumes that the data lies
 * in the z=0 plane.
 *
 *
*/

#ifndef vtkCookieCutFilter_h
#define vtkCookieCutFilter_h

#include "vtkFiltersModelingModule.h" // For export macro
#include "vtkPolyDataAlgorithm.h"

class VTKFILTERSMODELING_EXPORT vtkCookieCutFilter : public vtkPolyDataAlgorithm
{
public:
  //@{
  /**
   * Standard methods to instantiate, print and provide type information.
   */
  static vtkCookieCutFilter *New();
  vtkTypeMacro(vtkCookieCutFilter,vtkPolyDataAlgorithm);
  void PrintSelf(ostream& os, vtkIndent indent) VTK_OVERRIDE;
  //@}

  //@{
  /**
   * Set/Get the array of point coordinates defining the loop. There must
   * be at least three points used to define a loop.
   */
  virtual void SetLoop(vtkPoints*);
  vtkGetObjectMacro(Loop,vtkPoints);
  //@}

  /**
   * Overload GetMTime() because we depend on an input loop
   */
  vtkMTimeType GetMTime();

protected:
  vtkCookieCutFilter();
  ~vtkCookieCutFilter();

  int RequestData(vtkInformation *, vtkInformationVector **,
                  vtkInformationVector *) VTK_OVERRIDE;

  vtkPoints *Loop;

private:
  vtkCookieCutFilter(const vtkCookieCutFilter&) VTK_DELETE_FUNCTION;
  void operator=(const vtkCookieCutFilter&) VTK_DELETE_FUNCTION;
};


#endif
