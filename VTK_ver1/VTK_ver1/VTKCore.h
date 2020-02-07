#pragma once
#ifndef CORE_H
#define CORE_H

#include <vtkAutoInit.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2);
VTK_MODULE_INIT(vtkInteractionStyle);
VTK_MODULE_INIT(vtkRenderingFreeType);


//Sylinder Example
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkCylinderSource.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>

#include <array>


//Pointer Example
#include <vtkSmartPointer.h>

#include <vtkActor.h>
#include <vtkGlyph3D.h>
#include <vtkGlyphSource2D.h>
#include <vtkGraphLayout.h>
#include <vtkGraphLayoutView.h>
#include <vtkGraphToPolyData.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSimple2DLayoutStrategy.h>

//ReadSTL
#include <vtkBYUReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>
#include <vtksys/SystemTools.hxx>

//WriteSTL
#include <vtkPolyData.h>
#include <vtkSTLWriter.h>
#include <vtkSTLReader.h>
#include <vtkSphereSource.h>
#include <vtkSmartPointer.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

//Boolean
#include <vtkSphere.h>
#include <vtkBooleanOperationPolyDataFilter.h>
#include <vtkBox.h>
#include <vtkImplicitBoolean.h>
#include <vtkSampleFunction.h>
#include <vtkContourFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkNamedColors.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkProperty.h>
#include <vtkColor.h>

//BTEest

#include <vtkSmartPointer.h>
#include <vtkBooleanOperationPolyDataFilter.h>

#include <vtkCleanPolyData.h>
#include <vtkTriangleFilter.h>
#include <vtkActor.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkNamedColors.h>

#include <vtkBYUReader.h>
#include <vtkOBJReader.h>
#include <vtkPLYReader.h>
#include <vtkPolyDataReader.h>
#include <vtkSTLReader.h>
#include <vtkXMLPolyDataReader.h>

#include <vtkSphereSource.h>
#include <vtksys/SystemTools.hxx>

#include <vtkCamera.h>

//Subdivision
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkPointData.h>
#include <vtkSphereSource.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>

#include <string>

#include <vtkButterflySubdivisionFilter.h>
#include <vtkLoopSubdivisionFilter.h>
#include <vtkLinearSubdivisionFilter.h>

//Animation
#include "AnimateActors.h"

#include <vtkAnimationCue.h>
#include <vtkAnimationScene.h>
#include <vtkCamera.h>
#include <vtkColor.h>
#include <vtkCommand.h>
#include <vtkConeSource.h>
#include <vtkNamedColors.h>
#include <vtkPolyDataMapper.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>

class VTKCore
{
private:
	vtkSmartPointer<vtkSTLReader> Mesh =
		vtkSmartPointer<vtkSTLReader>::New();
protected:


public:
	VTKCore();
	~VTKCore();

	void SylinderExample();
	void PointerExample();
	void PositionCamera(vtkSmartPointer<vtkRenderer> &renderer,
		double *viewUp,
		double *position);
	vtkSmartPointer<vtkPolyData> ReadFile(std::string inputFilename);
	void ReadSTL(std::string inputFilename);
	void WriteSTL(std::string inputFilename);

	void DistanceBetweenPoints();

	void BooleanDemo();
	void BooleanTest();



	void AnimationExample();
};






#endif