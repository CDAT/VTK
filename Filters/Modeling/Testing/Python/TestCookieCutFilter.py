#!/usr/bin/env python

import vtk
from vtk.util.misc import vtkGetDataRoot

# create planes
# Create the RenderWindow, Renderer
#
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer( ren )

iren = vtk.vtkRenderWindowInteractor()
iren.SetRenderWindow(renWin)

# create pipeline. Use a plane source to create an array
# of points; then glyph them. Then cookie cut the glphs.
#
plane = vtk.vtkPlaneSource()
plane.SetXResolution(25)
plane.SetYResolution(25)

# Custom glyph
glyphData = vtk.vtkPolyData()
glyphPts = vtk.vtkPoints()
glyphVerts = vtk.vtkCellArray()
glyphLines = vtk.vtkCellArray()
glyphPolys = vtk.vtkCellArray()
glyphData.SetPoints(glyphPts)
glyphData.SetVerts(glyphVerts)
glyphData.SetLines(glyphLines)
glyphData.SetPolys(glyphPolys)

glyphPts.InsertPoint(0, -0.5,-0.5,0.0)
glyphPts.InsertPoint(1,  0.5,-0.5,0.0)
glyphPts.InsertPoint(2,  0.5, 0.5,0.0)
glyphPts.InsertPoint(3, -0.5, 0.5,0.0)
glyphPts.InsertPoint(4,  0.0,-0.5,0.0)
glyphPts.InsertPoint(5,  0.0,-0.75,0.0)
glyphPts.InsertPoint(6,  0.5, 0.0,0.0)
glyphPts.InsertPoint(7,  0.75, 0.0,0.0)
glyphPts.InsertPoint(8,  0.0, 0.5,0.0)
glyphPts.InsertPoint(9,  0.0, 0.75,0.0)
glyphPts.InsertPoint(10, -0.5,0.0,0.0)
glyphPts.InsertPoint(11, -0.75,0.0,0.0)
glyphPts.InsertPoint(12, 0.0,0.0,0.0)

glyphVerts.InsertNextCell(1)
glyphVerts.InsertCellPoint(12)

glyphLines.InsertNextCell(2)
glyphLines.InsertCellPoint(4)
glyphLines.InsertCellPoint(5)
glyphLines.InsertNextCell(2)
glyphLines.InsertCellPoint(6)
glyphLines.InsertCellPoint(7)
glyphLines.InsertNextCell(2)
glyphLines.InsertCellPoint(8)
glyphLines.InsertCellPoint(9)
glyphLines.InsertNextCell(2)
glyphLines.InsertCellPoint(10)
glyphLines.InsertCellPoint(11)

glyphPolys.InsertNextCell(4)
glyphPolys.InsertCellPoint(0)
glyphPolys.InsertCellPoint(1)
glyphPolys.InsertCellPoint(2)
glyphPolys.InsertCellPoint(3)

glyph = vtk.vtkGlyph3D()
glyph.SetInputConnection(plane.GetOutputPort())
glyph.SetSourceData(glyphData)
glyph.SetScaleFactor( 0.02 )

# Create a loop for cookie cutting
loop = vtk.vtkPoints()
loop.SetNumberOfPoints(4)
loop.SetPoint(0, -0.35,0.0,0.0)
loop.SetPoint(1, 0,-0.35,0.0)
loop.SetPoint(2, 0.35,0.0,0.0)
loop.SetPoint(3, 0.0,0.35,0.0)

cookie = vtk.vtkCookieCutFilter()
cookie.SetInputConnection(glyph.GetOutputPort())
cookie.SetLoop(loop)

mapper = vtk.vtkPolyDataMapper()
#mapper.SetInputConnection(glyph.GetOutputPort())
mapper.SetInputConnection(cookie.GetOutputPort())

actor = vtk.vtkActor()
actor.SetMapper(mapper)

loopData = vtk.vtkPolyData()
loopPoly = vtk.vtkCellArray()
loopPoly.InsertNextCell(4)
loopPoly.InsertCellPoint(0)
loopPoly.InsertCellPoint(1)
loopPoly.InsertCellPoint(2)
loopPoly.InsertCellPoint(3)
loopData.SetPoints(loop)
loopData.SetPolys(loopPoly)

loopMapper = vtk.vtkPolyDataMapper()
loopMapper.SetInputData(loopData)

loopActor = vtk.vtkActor()
loopActor.SetMapper(loopMapper)
loopActor.GetProperty().SetColor(1,0,0)
loopActor.GetProperty().SetRepresentationToWireframe()

ren.AddActor(actor)
ren.AddActor(loopActor)

renWin.Render()
iren.Start()
