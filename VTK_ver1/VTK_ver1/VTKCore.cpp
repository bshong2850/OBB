#include "VTKCore.h"

VTKCore::VTKCore()
{

}

VTKCore::~VTKCore()
{

}

void VTKCore::SylinderExample()
{
	vtkSmartPointer<vtkNamedColors> colors =
		vtkSmartPointer<vtkNamedColors>::New();
		
	// Set the background color.
	std::array<unsigned char, 4> bkg{ { 26, 51, 102, 255 } };
	colors->SetColor("BkgColor", bkg.data());
		
	// This creates a polygonal cylinder model with eight circumferential facets
	// (i.e, in practice an octagonal prism).
	vtkSmartPointer<vtkCylinderSource> cylinder =
		vtkSmartPointer<vtkCylinderSource>::New();
	cylinder->SetResolution(8);
		
	// The mapper is responsible for pushing the geometry into the graphics
	// library. It may also do color mapping, if scalars or other attributes are
	// defined.
	vtkSmartPointer<vtkPolyDataMapper> cylinderMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	cylinderMapper->SetInputConnection(cylinder->GetOutputPort());
		
	// The actor is a grouping mechanism: besides the geometry (mapper), it
	// also has a property, transformation matrix, and/or texture map.
	// Here we set its color and rotate it around the X and Y axes.
	vtkSmartPointer<vtkActor> cylinderActor = vtkSmartPointer<vtkActor>::New();
	cylinderActor->SetMapper(cylinderMapper);
	cylinderActor->GetProperty()->SetColor(
		colors->GetColor4d("Tomato").GetData());
	cylinderActor->RotateX(30.0);
	cylinderActor->RotateY(-45.0);
		
	// The renderer generates the image
	// which is then displayed on the render window.
	// It can be thought of as a scene to which the actor is added
	vtkSmartPointer<vtkRenderer> renderer = vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(cylinderActor);
	renderer->SetBackground(colors->GetColor3d("BkgColor").GetData());
	// Zoom in a little by accessing the camera and invoking its "Zoom" method.
	renderer->ResetCamera();
	renderer->GetActiveCamera()->Zoom(1.5);
		
	// The render window is the actual GUI window
	// that appears on the computer screen
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(600, 600);
	renderWindow->AddRenderer(renderer);
	renderWindow->SetWindowName("Cylinder");
		
	// The render window interactor captures mouse events
	// and will perform appropriate camera or actor manipulation
	// depending on the nature of the events.
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
		
	// This starts the event loop and as a side effect causes an initial render.
	renderWindow->Render();
	renderWindowInteractor->Start();
		
	return ;
}

void VTKCore::PointerExample()
{
	vtkSmartPointer<vtkMutableDirectedGraph> g =
		vtkSmartPointer<vtkMutableDirectedGraph>::New();

	vtkIdType v1 = g->AddVertex();
	vtkIdType v2 = g->AddVertex();
	vtkIdType v3 = g->AddVertex();

	g->AddEdge(v1, v2);
	g->AddEdge(v2, v3);
	g->AddEdge(v3, v1);

	// Do layout manually before handing graph to the view.
	// This allows us to know the positions of edge arrows.
	vtkSmartPointer<vtkGraphLayoutView> graphLayoutView =
		vtkSmartPointer<vtkGraphLayoutView>::New();

	vtkSmartPointer<vtkGraphLayout> layout =
		vtkSmartPointer<vtkGraphLayout>::New();
	vtkSmartPointer<vtkSimple2DLayoutStrategy> strategy =
		vtkSmartPointer<vtkSimple2DLayoutStrategy>::New();
	layout->SetInputData(g);
	layout->SetLayoutStrategy(strategy);

	// Tell the view to use the vertex layout we provide
	graphLayoutView->SetLayoutStrategyToPassThrough();
	// The arrows will be positioned on a straight line between two
	// vertices so tell the view not to draw arcs for parallel edges
	graphLayoutView->SetEdgeLayoutStrategyToPassThrough();

	// Add the graph to the view. This will render vertices and edges,
	// but not edge arrows.
	graphLayoutView->AddRepresentationFromInputConnection(
		layout->GetOutputPort());

	// Manually create an actor containing the glyphed arrows.
	vtkSmartPointer<vtkGraphToPolyData> graphToPoly = vtkSmartPointer<vtkGraphToPolyData>::New();
	graphToPoly->SetInputConnection(layout->GetOutputPort());
	graphToPoly->EdgeGlyphOutputOn();

	// Set the position (0: edge start, 1: edge end) where
	// the edge arrows should go.
	graphToPoly->SetEdgeGlyphPosition(0.98);

	// Make a simple edge arrow for glyphing.
	vtkSmartPointer<vtkGlyphSource2D> arrowSource =	vtkSmartPointer<vtkGlyphSource2D>::New();
	arrowSource->SetGlyphTypeToEdgeArrow();
	arrowSource->SetScale(0.1);
	arrowSource->Update();

	// Use Glyph3D to repeat the glyph on all edges.
	vtkSmartPointer<vtkGlyph3D> arrowGlyph =
		vtkSmartPointer<vtkGlyph3D>::New();
	arrowGlyph->SetInputConnection(0, graphToPoly->GetOutputPort(1));
	arrowGlyph->SetInputConnection(1, arrowSource->GetOutputPort());

	// Add the edge arrow actor to the view.
	vtkSmartPointer<vtkPolyDataMapper> arrowMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	arrowMapper->SetInputConnection(arrowGlyph->GetOutputPort());
	vtkSmartPointer<vtkActor> arrowActor =
		vtkSmartPointer<vtkActor>::New();
	arrowActor->SetMapper(arrowMapper);
	graphLayoutView->GetRenderer()->AddActor(arrowActor);

	graphLayoutView->ResetCamera();
	graphLayoutView->Render();
	graphLayoutView->GetInteractor()->Start();

	return ;
}
void VTKCore::PositionCamera(vtkSmartPointer<vtkRenderer> &renderer,
	double *viewUp,
	double *position)
{
	renderer->GetActiveCamera()->SetFocalPoint(0.0, 0.0, 0.0);
	renderer->GetActiveCamera()->SetViewUp(viewUp);
	renderer->GetActiveCamera()->SetPosition(position);
	renderer->ResetCamera();
	return;
}
void VTKCore::ReadSTL(std::string inputFilename)
{
	
	Mesh->SetFileName(inputFilename.c_str());
	Mesh->Update();

	// Visualize
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(Mesh->GetOutputPort());

	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	vtkSmartPointer<vtkRenderWindowInteractor> renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);

	renderer->AddActor(actor);
	renderer->SetBackground(.3, .6, .3); // Background color green

	renderWindow->Render();
	renderWindowInteractor->Start();

	return ;
}

vtkSmartPointer<vtkPolyData> VTKCore::ReadFile(std::string inputFilename)
{
	vtkSmartPointer<vtkPolyData> polyData;
	std::string extension = vtksys::SystemTools::GetFilenameExtension(inputFilename);
	if (extension == ".ply")
	{
		vtkSmartPointer<vtkPLYReader> reader =
			vtkSmartPointer<vtkPLYReader>::New();
		reader->SetFileName(inputFilename.c_str());
		reader->Update();
		polyData = reader->GetOutput();
	}
	else if (extension == ".vtp")
	{
		vtkSmartPointer<vtkXMLPolyDataReader> reader =
			vtkSmartPointer<vtkXMLPolyDataReader>::New();
		reader->SetFileName(inputFilename.c_str());
		reader->Update();
		polyData = reader->GetOutput();
	}
	else if (extension == ".obj")
	{
		vtkSmartPointer<vtkOBJReader> reader =
			vtkSmartPointer<vtkOBJReader>::New();
		reader->SetFileName(inputFilename.c_str());
		reader->Update();
		polyData = reader->GetOutput();
	}
	else if (extension == ".stl")
	{
		vtkSmartPointer<vtkSTLReader> reader =
			vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(inputFilename.c_str());
		reader->Update();
		polyData = reader->GetOutput();
	}
	else if (extension == ".vtk")
	{
		vtkSmartPointer<vtkPolyDataReader> reader =
			vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(inputFilename.c_str());
		reader->Update();
		polyData = reader->GetOutput();
	}
	else if (extension == ".g")
	{
		vtkSmartPointer<vtkBYUReader> reader =
			vtkSmartPointer<vtkBYUReader>::New();
		reader->SetGeometryFileName(inputFilename.c_str());
		reader->Update();
		polyData = reader->GetOutput();
	}
	else
	{
		vtkSmartPointer<vtkSphereSource> source =
			vtkSmartPointer<vtkSphereSource>::New();
		source->Update();
		polyData = source->GetOutput();
	}
	return polyData;
}

void VTKCore::WriteSTL(std::string writeFilename)
{
	vtkSmartPointer<vtkSTLWriter> stlWriter =
		vtkSmartPointer<vtkSTLWriter>::New();
	stlWriter->SetFileName(writeFilename.c_str());
	stlWriter->SetInputConnection(Mesh->GetOutputPort());
	stlWriter->Write();
}

void VTKCore::DistanceBetweenPoints()
{
	double p0[3] = { 0.0, 0.0, 0.0 };
	double p1[3] = { 1.0, 1.0, 1.0 };

	// Find the squared distance between the points.
	double squaredDistance = vtkMath::Distance2BetweenPoints(p0, p1);

	// Take the square root to get the Euclidean distance between the points.
	double distance = sqrt(squaredDistance);

	// Output the results.
	std::cout << "SquaredDistance = " << squaredDistance << std::endl;
	std::cout << "Distance = " << distance << std::endl;

	return ;
}

void VTKCore::BooleanDemo()
{
	// Define colors
	auto colors =
		vtkSmartPointer<vtkNamedColors>::New();
	vtkColor3d actorColor = colors->GetColor3d("AliceBlue");
	vtkColor3d  EdgeColour = colors->GetColor3d("SteelBlue");
	vtkColor3d BackgroundColour = colors->GetColor3d("Silver");

	std::string inputFilename = "model/basic_cube.stl";
	Mesh->SetFileName(inputFilename.c_str());
	Mesh->Update();
	auto Mesh_ = vtkSmartPointer<vtkPolyDataMapper>::New();
	Mesh_->SetInputConnection(Mesh->GetOutputPort());

	// create a sphere
	auto sphere =
		vtkSmartPointer<vtkSphere>::New();
	sphere->SetCenter(0.0, 0.0, 0.0);
	sphere->SetRadius(1);

	//create a box
	auto box =
		vtkSmartPointer<vtkBox>::New();
	box->SetBounds(-1, 1, -1, 1, -1, 1);

	// combine the two implicit functions
	auto boolean =
		vtkSmartPointer<vtkImplicitBoolean>::New();
	boolean->SetOperationTypeToDifference();
	/*boolean->SetInputConnection(0, Mesh->GetOutputPort());
	boolean->SetInputConnection(1, sphere->GetOutputPort());*/
	boolean->AddFunction(box);
	boolean->AddFunction(sphere);
	//boolean->SetOperationTypeToDifference();
	/*auto boolean =
		vtkSmartPointer<vtkImplicitBoolean>::New();
	boolean->SetOperationTypeToDifference();*/
	// boolean->SetOperationTypeToUnion()
	// boolean->SetOperationTypeToIntersection()
	/*boolean->AddFunction(box);
	boolean->AddFunction(sphere);*/
	//The sample function generates a distance function from the implicit
	//function.This is then contoured to get a polygonal surface.
	auto sample =
		vtkSmartPointer<vtkSampleFunction>::New();
	sample->SetImplicitFunction(boolean);
	sample->SetModelBounds(-1, 2, -1, 1, -1, 1);
	sample->SetSampleDimensions(20, 20, 20);
	sample->ComputeNormalsOff();

	// contour
	auto surface =
		vtkSmartPointer<vtkContourFilter>::New();
	surface->SetInputConnection(sample->GetOutputPort());
	surface->SetValue(0, 0.0);

	//Create a mapper and an actor
	auto mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(surface->GetOutputPort());
	mapper->ScalarVisibilityOff();
	auto actor =
		vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->EdgeVisibilityOn();
	actor->GetProperty()->SetColor(actorColor.GetData());
	actor->GetProperty()->SetEdgeColor(EdgeColour.GetData());

	// A renderer and render window
	auto renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->SetBackground(BackgroundColour.GetData());
	auto renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->SetSize(600, 600);
	renderWindow->AddRenderer(renderer);
	auto renderWindowInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renderWindowInteractor->SetRenderWindow(renderWindow);
	//add the actor
	renderer->AddActor(actor);
	// Start
	renderer->GetActiveCamera()->SetPosition(5.0, -4.0, 1.6);
	renderer->GetActiveCamera()->SetViewUp(0.1, 0.5, 0.9);
	renderer->GetActiveCamera()->SetDistance(6.7);
	renderWindow->Render();
	renderWindowInteractor->Start();

	return;
}

void VTKCore::BooleanTest()
{
	vtkSmartPointer<vtkPolyData> input1;
	vtkSmartPointer<vtkPolyData> input2;

	std::string operation("difference");
	//union
	//difference
	//intersection
	//std::string operation("intersection");
	vtkSmartPointer<vtkPolyData> poly1;
	//poly1 = ReadFile("model/armadillo.obj");
	poly1 = ReadFile("model/basic_cube_Sub_3.stl");
	vtkSmartPointer<vtkTriangleFilter> tri1 =
		vtkSmartPointer<vtkTriangleFilter>::New();
	tri1->SetInputData(poly1);
	vtkSmartPointer<vtkCleanPolyData> clean1 =
		vtkSmartPointer<vtkCleanPolyData>::New();
	clean1->SetInputConnection(tri1->GetOutputPort());
	clean1->Update();
	input1 = clean1->GetOutput();

	/*vtkSmartPointer<vtkPolyDataAlgorithm> subdivisionFilter_1;
	subdivisionFilter_1 = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
	dynamic_cast<vtkLinearSubdivisionFilter *> (subdivisionFilter_1.GetPointer())->SetNumberOfSubdivisions(5);
	subdivisionFilter_1->SetInputData(input1);
	subdivisionFilter_1->Update();*/

	vtkSmartPointer<vtkPolyData> poly2;
	poly2 = ReadFile("model/basic_cube_Sub_3.stl");
	vtkSmartPointer<vtkTriangleFilter> tri2 =
		vtkSmartPointer<vtkTriangleFilter>::New();
	tri2->SetInputData(poly2);
	tri2->Update();
	vtkSmartPointer<vtkCleanPolyData> clean2 =
		vtkSmartPointer<vtkCleanPolyData>::New();
	clean2->SetInputConnection(tri2->GetOutputPort());
	clean2->Update();
	input2 = clean2->GetOutput();

	/*vtkSmartPointer<vtkPolyDataAlgorithm> subdivisionFilter_2;
	subdivisionFilter_2 = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
	dynamic_cast<vtkLinearSubdivisionFilter *> (subdivisionFilter_2.GetPointer())->SetNumberOfSubdivisions(5);
	subdivisionFilter_2->SetInputData(input2);
	subdivisionFilter_2->Update();*/



	/*switch (i)
	{
	case 0:
	subdivisionFilter = vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
	dynamic_cast<vtkLinearSubdivisionFilter *> (subdivisionFilter.GetPointer())->SetNumberOfSubdivisions(3);
	break;
	case 1:
	subdivisionFilter = vtkSmartPointer<vtkLoopSubdivisionFilter>::New();
	dynamic_cast<vtkLoopSubdivisionFilter *> (subdivisionFilter.GetPointer())->SetNumberOfSubdivisions(3);
	break;
	case 2:
	subdivisionFilter = vtkSmartPointer<vtkButterflySubdivisionFilter>::New();
	dynamic_cast<vtkButterflySubdivisionFilter *> (subdivisionFilter.GetPointer())->SetNumberOfSubdivisions(3);
	break;
	default:
	break;
	}
	*/

	/*vtkSmartPointer<vtkTriangleFilter> triangles =
	vtkSmartPointer<vtkTriangleFilter>::New();
	triangles->SetInputConnection(tri2->GetOutputPort());
	triangles->Update();
	input2 = triangles->GetOutput();*/
	/*auto sample =
	vtkSmartPointer<vtkSampleFunction>::New();
	sample->SetImplicitFunction(input1);
	sample->SetModelBounds(-1, 2, -1, 1, -1, 1);
	sample->SetSampleDimensions(40, 40, 40);
	sample->ComputeNormalsOff();*/

	//operation = 2;
	/* if (argc == 4)
	{
	vtkSmartPointer<vtkPolyData> poly1;
	poly1 = ReadPolyData(argv[1]);
	vtkSmartPointer<vtkTriangleFilter> tri1 =
	vtkSmartPointer<vtkTriangleFilter>::New();
	tri1->SetInputData(poly1);
	vtkSmartPointer<vtkCleanPolyData> clean1 =
	vtkSmartPointer<vtkCleanPolyData>::New();
	clean1->SetInputConnection(tri1->GetOutputPort());
	clean1->Update();
	input1 = clean1->GetOutput();

	vtkSmartPointer<vtkPolyData> poly2;
	poly2 = ReadPolyData(argv[3]);
	vtkSmartPointer<vtkTriangleFilter> tri2 =
	vtkSmartPointer<vtkTriangleFilter>::New();
	tri2->SetInputData(poly2);
	tri2->Update();
	vtkSmartPointer<vtkCleanPolyData> clean2 =
	vtkSmartPointer<vtkCleanPolyData>::New();
	clean2->SetInputConnection(tri2->GetOutputPort());
	clean2->Update();
	input2 = clean2->GetOutput();
	operation = argv[2];
	}*/
	/* else
	{
	vtkSmartPointer<vtkSphereSource> sphereSource1 =
	vtkSmartPointer<vtkSphereSource>::New();
	sphereSource1->SetCenter(.25, 0, 0);
	sphereSource1->SetPhiResolution(21);
	sphereSource1->SetThetaResolution(21);
	sphereSource1->Update();
	input1 = sphereSource1->GetOutput();

	vtkSmartPointer<vtkSphereSource> sphereSource2 =
	vtkSmartPointer<vtkSphereSource>::New();
	sphereSource2->Update();
	input2 = sphereSource2->GetOutput();

	if (argc == 2)
	{
	operation = argv[1];
	}
	}*/

	vtkSmartPointer<vtkNamedColors> colors =
		vtkSmartPointer<vtkNamedColors>::New();

	vtkSmartPointer<vtkPolyDataMapper> input1Mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	input1Mapper->SetInputData(input1);
	//input1Mapper->SetInputConnection(subdivisionFilter_1->GetOutputPort());
	input1Mapper->ScalarVisibilityOff();
	vtkSmartPointer<vtkActor> input1Actor =
		vtkSmartPointer<vtkActor>::New();
	input1Actor->SetMapper(input1Mapper);
	input1Actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Tomato").GetData());
	input1Actor->GetProperty()->SetSpecular(.6);
	input1Actor->GetProperty()->SetSpecularPower(5);


	input1Actor->SetPosition(input1->GetBounds()[1] - input1->GetBounds()[0], 0, 0);
	vtkSmartPointer<vtkPolyDataMapper> input2Mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	input2Mapper->SetInputData(input2);
	//input2Mapper->SetInputConnection(subdivisionFilter_2->GetOutputPort());
	input2Mapper->ScalarVisibilityOff();
	vtkSmartPointer<vtkActor> input2Actor =
		vtkSmartPointer<vtkActor>::New();
	input2Actor->SetMapper(input2Mapper);
	input2Actor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Mint").GetData());
	input2Actor->GetProperty()->SetSpecular(.6);
	input2Actor->GetProperty()->SetSpecularPower(5);
	input2Actor->SetPosition(
		-(input2->GetBounds()[1] - input2->GetBounds()[0]),
		0, 0);


	vtkColor3d  EdgeColour = colors->GetColor3d("SteelBlue");
	/*input1Actor->GetProperty()->EdgeVisibilityOn();
	input1Actor->GetProperty()->SetEdgeColor(EdgeColour.GetData());*/
	//vtkColor3d  EdgeColour = colors->GetColor3d("SteelBlue");
	input2Actor->GetProperty()->EdgeVisibilityOn();
	input2Actor->GetProperty()->SetEdgeColor(EdgeColour.GetData());


	vtkSmartPointer<vtkBooleanOperationPolyDataFilter> booleanOperation =
		vtkSmartPointer<vtkBooleanOperationPolyDataFilter>::New();
	if (operation == "union")
	{
		booleanOperation->SetOperationToUnion();
	}
	else if (operation == "intersection")
	{
		booleanOperation->SetOperationToIntersection();
	}
	else if (operation == "difference")
	{
		booleanOperation->SetOperationToDifference();
	}
	else
	{
		std::cout << "Unknown operation: " << operation << std::endl;
		return ;
	}
	booleanOperation->SetInputData(0, input1);
	booleanOperation->SetInputData(1, input2);
	/*booleanOperation->SetInputConnection(0, subdivisionFilter_1->GetOutputPort());
	booleanOperation->SetInputConnection(1, subdivisionFilter_2->GetOutputPort());*/

	vtkSmartPointer<vtkPolyDataMapper> booleanOperationMapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	booleanOperationMapper->SetInputConnection(booleanOperation->GetOutputPort());
	booleanOperationMapper->ScalarVisibilityOff();

	vtkSmartPointer<vtkActor> booleanOperationActor =
		vtkSmartPointer<vtkActor>::New();
	booleanOperationActor->SetMapper(booleanOperationMapper);
	booleanOperationActor->GetProperty()->SetDiffuseColor(colors->GetColor3d("Banana").GetData());
	booleanOperationActor->GetProperty()->SetSpecular(.6);
	booleanOperationActor->GetProperty()->SetSpecularPower(5);



	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->AddViewProp(input1Actor);
	renderer->AddViewProp(input2Actor);
	renderer->AddViewProp(booleanOperationActor);
	renderer->SetBackground(colors->GetColor3d("Silver").GetData());
	vtkSmartPointer<vtkRenderWindow> renderWindow =
		vtkSmartPointer<vtkRenderWindow>::New();
	renderWindow->AddRenderer(renderer);
	renderWindow->SetSize(640, 480);
	double viewUp[3] = { 0.0, 0.0, 1.0 };
	double position[3] = { 0.0, -1.0, 0.0 };
	PositionCamera(renderer, viewUp, position);
	renderer->GetActiveCamera()->Dolly(1.5);
	renderer->ResetCameraClippingRange();

	vtkSmartPointer<vtkRenderWindowInteractor> renWinInteractor =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	renWinInteractor->SetRenderWindow(renderWindow);

	renderWindow->Render();
	renWinInteractor->Start();

	return ;
}

void VTKCore::AnimationExample()
{
	// Colors
	auto colors = vtkSmartPointer<vtkNamedColors>::New();
	vtkColor3d coneColor = colors->GetColor3d("Tomato");
	vtkColor3d sphereColor = colors->GetColor3d("Banana");
	vtkColor3d backgroundColor = colors->GetColor3d("Peacock");

	// Create the graphics structure. The renderer renders into the
	// render window.
	auto iren = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	auto ren1 = vtkSmartPointer<vtkRenderer>::New();
	ren1->SetBackground(backgroundColor.GetData());

	auto renWin = vtkSmartPointer<vtkRenderWindow>::New();
	renWin->SetMultiSamples(0);
	iren->SetRenderWindow(renWin);
	renWin->AddRenderer(ren1);

	// Generate a sphere
	auto sphereSource = vtkSmartPointer<vtkSphereSource>::New();
	sphereSource->SetPhiResolution(31);
	sphereSource->SetThetaResolution(31);

	auto sphereMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	sphereMapper->SetInputConnection(sphereSource->GetOutputPort());
	auto sphere = vtkSmartPointer<vtkActor>::New();
	sphere->SetMapper(sphereMapper);
	sphere->GetProperty()->SetDiffuseColor(sphereColor.GetData());
	sphere->GetProperty()->SetDiffuse(.7);
	sphere->GetProperty()->SetSpecular(.3);
	sphere->GetProperty()->SetSpecularPower(30.0);

	ren1->AddActor(sphere);

	// Generate a cone
	auto coneSource = vtkSmartPointer<vtkConeSource>::New();
	coneSource->SetResolution(31);

	auto coneMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	coneMapper->SetInputConnection(coneSource->GetOutputPort());
	auto cone = vtkSmartPointer<vtkActor>::New();
	cone->SetMapper(coneMapper);
	cone->GetProperty()->SetDiffuseColor(coneColor.GetData());

	ren1->AddActor(cone);

	// Create an Animation Scene
	auto scene = vtkSmartPointer<vtkAnimationScene>::New();
	/*if (argc >= 2 && strcmp(argv[1], "-real") == 0)
	{
		cout << "real-time mode" << endl;
		scene->SetModeToRealTime();
	}
	else
	{
		cout << "sequence mode" << endl;
		
	}*/
	scene->SetModeToRealTime();
	//scene->SetModeToSequence();
	scene->SetLoop(0);
	scene->SetFrameRate(5);
	scene->SetStartTime(0);
	scene->SetEndTime(20);

	vtkSmartPointer<AnimationSceneObserver> sceneObserver =
		vtkSmartPointer<AnimationSceneObserver>::New();
	sceneObserver->SetRenderWindow(renWin);
	scene->AddObserver(vtkCommand::AnimationCueTickEvent, sceneObserver);

	// Create an Animation Cue for each actor
	auto cue1 = vtkSmartPointer<vtkAnimationCue>::New();
	cue1->SetStartTime(5);
	cue1->SetEndTime(23);
	scene->AddCue(cue1);

	auto cue2 = vtkSmartPointer<vtkAnimationCue>::New();
	cue2->SetStartTime(1);
	cue2->SetEndTime(10);
	scene->AddCue(cue2);

	// Create an ActorAnimator for each actor;
	ActorAnimator animateSphere;
	animateSphere.SetActor(sphere);
	animateSphere.AddObserversToCue(cue1);

	ActorAnimator animateCone;
	std::vector<double> endCone(3);
	endCone[0] = -1;
	endCone[1] = -1;
	endCone[2] = -1;
	animateCone.SetEndPosition(endCone);
	animateCone.SetActor(cone);
	animateCone.AddObserversToCue(cue2);

	renWin->Render();
	ren1->ResetCamera();
	ren1->GetActiveCamera()->Dolly(.5);
	ren1->ResetCameraClippingRange();

	// Create Cue observer.
	scene->Play();
	scene->Stop();

	iren->Start();
	return ;
}