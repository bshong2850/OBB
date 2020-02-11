#pragma once

#include "iglcore.h"

IGLCore::IGLCore()
{
	Initialize();
}

IGLCore::~IGLCore()
{

}

void IGLCore::Initialize()
{
	std::string full_path = p_writefilepath;

	char temp[256];
	strcpy(temp, full_path.c_str());
	char *p = temp;
	while (*p)
	{
		if (('\\' == *p) || ('/' == *p))
		{
			if (':' != *(p - 1))
			{
				CreateDirectory(temp, NULL);
			}
		}
		*p++;
	}
	_v = Eigen::MatrixXd::Zero(1, 1);
	_f = Eigen::MatrixXi::Zero(1, 1);
	_n = Eigen::MatrixXd::Zero(1, 1);
}

void IGLCore::PrintInformation()
{
	if(_v.size() > 3)
	{
		std::cout << "Vertex Information" << std::endl;
		std::cout << " Vertex Size = " << _v.cols() << " " << _v.rows() << std::endl;
		for (int i = 0; i < _v.rows(); i++)
		{
			std::cout << "   ";
			for (int j = 0; j < _v.cols(); j++)
			{
				std::cout.width(10);
				std::cout << std::left << _v(i, j) << " ";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "Vertex is EMPTY" << std::endl;
	}
	if (_f.size() > 1)
	{
		std::cout << "Face Information" << std::endl;
		std::cout << " Face Size = " << _f.cols() << " " << _f.rows() << std::endl;
		for (int i = 0; i < _f.rows(); i++)
		{
			std::cout << "   ";
			for (int j = 0; j < _f.cols(); j++)
			{
				std::cout.width(2);
				std::cout << std::left << _f(i, j) << " ";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "Face is EMPTY" << std::endl;
	}
	if (_n.size() > 1)
	{
		std::cout << "Normal Information" << std::endl;
		std::cout << " Normal Size = " << _n.cols() << " " << _n.rows() << std::endl;
		for (int i = 0; i < _n.rows(); i++)
		{
			std::cout << "   ";
			for (int j = 0; j < _n.cols(); j++)
			{
				std::cout.width(10);
				std::cout << std::left << _n(i, j) << " ";
			}
			std::cout << std::endl;
		}
	}
	else
	{
		std::cout << "Normal is EMPTY" << std::endl;
	}
}

void IGLCore::ReadFile(std::string filename, std::string extension)
{
	std::string fullfilename = filename + "." + extension;
	if (extension == "stl")
		ReadSTL(fullfilename);
	else if (extension == "obj")
		ReadOBJ(fullfilename);
	else
	{
		std::cout << "This fileformat is not supported" << std::endl;
		exit(0);
	}
}
void IGLCore::ReadSTL(std::string filename, bool _duplicate)
{
	igl::readSTL(p_readfilepath + filename, p_originV, p_originF, p_originN);
	if (_duplicate)
	{
		Eigen::MatrixXd SI;
		Eigen::MatrixXd SVJ;
		igl::remove_duplicate_vertices(p_originV, 0.000001, _v, SI, SVJ);
		_f.resizeLike(p_originF);
		for (int f = 0;f < p_originF.rows();f++)
		{
			for (int c = 0;c < p_originF.cols();c++)
			{
				_f(f, c) = SVJ(p_originF(f, c));
			}
		}
		igl::per_face_normals_stable(_v, _f, _n);

		//Eigen::VectorXi _;
		//igl::bfs_orient(Eigen::MatrixXi(_f), _f, _);
		//igl::writeSTL("D:/hbs/Code/LibIGL_ver2/LibIGL_ver2/model/testestestsets.stl", _v, _f, _n, true);
	}
	else
	{
		_v = p_originV;
		_f = p_originF;
		_n = p_originN;
	}
	if (_v.size() >= 3)
		std::cout << "Read STL COMPLETE" << std::endl;
	else
		std::cout << "Read STL FAIL -> Mesh is Empty" << std::endl;
}

void IGLCore::ReadOBJ(std::string filename)
{
	igl::readOBJ(p_readfilepath + filename, p_originV, p_originF);
	if (p_originF.cols() >= 4)
	{
		igl::polygon_mesh_to_triangle_mesh(p_originF, _f);
		_v = p_originV;
	}
	else
	{
		_v = p_originV;
		_f = p_originF;
	}
	if (_v.size() >= 3)
		std::cout << "Read OBJ COMPLETE" << std::endl;
	else
		std::cout << "Read OBJ FAIL -> Mesh is Empty" << std::endl;
}

void IGLCore::WriteSTL(std::string filename, bool _ascii)
{
	if (_v.size() < 3)
	{
		std::cout << "Write STL FAIL -> Mesh is Empty" << std::endl;
		exit(0);
	}
	else
	{
		igl::writeSTL(p_writefilepath + filename, _v, _f, _n, _ascii);
		std::cout << "Write STL COMPLETE" << std::endl;
	}

}

void IGLCore::WriteDMAT(std::string filename)
{
	writeFile.open(p_readfilepath + filename);
	int _size = _v.size() / 3;
	int minus_count = 0;
	if (writeFile.is_open())
	{
		writeFile << "1" << " " << _size << std::endl;
		for (int i = 0; i < _size; i++)
		{
			//writeFile << "1" << std::endl;
			if (minus_count % 15 == 0)
			{
				writeFile << "-1" << std::endl;
			}
			else
			{
				writeFile << "1" << std::endl;
			}
			minus_count++;
		}
	}
	writeFile.close();
}

void IGLCore::ComputeNormal(int mode)
{
	if(VERTEX_NORMAL)
	{
		igl::per_vertex_normals(_v, _f, _n);
	}
	else if (FACE_NORMAL)
	{
		igl::per_face_normals(_v, _f, _n);
	}
	else if (CORNER_NORMAL)
	{
		igl::per_corner_normals(_v, _f, 20, _n);
	}
	std::cout << "Compute Normal COMPLETE" << std::endl;
}

void IGLCore::Smoothing()
{
	Eigen::SparseMatrix<double> L, G, K, M;
	Eigen::MatrixXd BC;
	Eigen::RowVector3d centroid(0, 0, 0);
	Eigen::VectorXd dblA;
	double area = 0;
	// Compute Laplace-Beltrami operator: #V by #V
	igl::cotmatrix(_v, _f, L);
	// Alternative construction of same Laplacian
	// Gradient/Divergence
	igl::grad(_v, _f, G);
	// Diagonal per-triangle "mass matrix"
	igl::doublearea(_v, _f, dblA);
	// Place areas along diagonal #dim times
	igl::massmatrix(_v, _f, igl::MASSMATRIX_TYPE_BARYCENTRIC, M);
	// Solve (M-delta*L) U = M*U
	const auto & S = (M - 0.001*L);
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver(S);
	assert(solver.info() == Eigen::Success);
	_v = solver.solve(M*_v).eval();
	// Compute centroid and subtract (also important for numerics)
	igl::doublearea(_v, _f, dblA);
	area = 0.5*dblA.sum();
	igl::barycenter(_v, _f, BC);
	for (int i = 0;i < BC.rows();i++)
	{
		centroid += 0.5*dblA(i) / area*BC.row(i);
	}
	_v.rowwise() -= centroid;
	_v.array() /= sqrt(area);
	//_v = U;
	//_v = U;
	std::cout << "Smoothing COMPLETE" << std::endl;
}

void IGLCore::Decimation()
{
	/*using namespace std;
	using namespace Eigen;
	using namespace igl;

	MatrixXd V, OV;
	MatrixXi F;
	VectorXi OF;
	igl::decimate(_v, _f, 100, V, F, OF);
	std::cout << "complete" << std::endl;*/
	using namespace std;
	using namespace Eigen;
	using namespace igl;
	// Load a closed manifold mesh
	
	MatrixXd V, OV;
	MatrixXi F, OF;
	
	V = _v;
	F = _f;
	// Prepare array-based edge data structures and priority queue
	VectorXi EMAP;
	MatrixXi E, EF, EI;
	typedef std::set<std::pair<double, int> > PriorityQueue;
	PriorityQueue Q;
	std::vector<PriorityQueue::iterator > Qit;
	// If an edge were collapsed, we'd collapse it to these points:
	MatrixXd C;
	int num_collapsed;

	// Function to reset original mesh and data structures
	edge_flaps(F, E, EMAP, EF, EI);
	Qit.resize(E.rows());

	C.resize(E.rows(), V.cols());
	VectorXd costs(E.rows());
	Q.clear();
	for (int e = 0;e < E.rows();e++)
	{
		double cost = e;
		RowVectorXd p(1, 3);
		shortest_edge_and_midpoint(e, V, F, E, EMAP, EF, EI, cost, p);
		C.row(e) = p;
		Qit[e] = Q.insert(std::pair<double, int>(cost, e)).first;
	}
	num_collapsed = 0;
	for (int i =0; i<10; i++)
	{
		bool something_collapsed = false;
		//const int max_iter = std::ceil(0.01*Q.size());
		const int max_iter = 10;
		for (int j = 0;j < max_iter;j++)
		{
			if (!collapse_edge(
				shortest_edge_and_midpoint, V, F, E, EMAP, EF, EI, Q, Qit, C))
			{
				break;
			}
			else
			{
				something_collapsed = true;
				Eigen::Matrix<int, Eigen::Dynamic, 1> I;
				/*igl::remove_duplicates(V, F, _v, _f, I);
				igl::per_face_normals(_v, _f, _n);
				V = _v;
				F = _f;*/
				num_collapsed++;
				break;
			}
		}
		if (something_collapsed)
		{
			_v = V;
			_f = F;
			WriteSTL("test_" + std::to_string(num_collapsed) + ".stl", true);
		}
		//bool something_collapsed = false;
		//// collapse edge
		//const int max_iter = std::ceil(0.01*Q.size());
		////const int max_iter = 15000;
		//for (int j = 0;j < max_iter;j++)
		//{
		//	//if (j % 100 == 0)
		//	std::cout << j << std::endl;
		//	if (collapse_edge(
		//		shortest_edge_and_midpoint, V, F, E, EMAP, EF, EI, Q, Qit, C))
		//	{
		//		if (num_collapsed % 10 == 0)
		//		{
		//			std::cout << "collapse!!!!" << std::endl;
		//			_v.resize(V.rows(), V.cols());
		//			_f.resize(F.rows(), F.cols());
		//			_v = V;
		//			_f = F;
		//			//igl::per_face_normals(_v, _f, _n);
		//			WriteSTL("Sub_3" + std::to_string(num_collapsed) + "_.stl");
		//		}
		//		num_collapsed++;
		//	}
		//}
		std::cout << "Decimation COMPLETE" << std::endl;
	}

}

void IGLCore::Subdivision()
{
	using namespace Eigen;
	using namespace igl;	

	igl::loop(_v, _f, _v, _f);
	igl::per_face_normals(_v, _f, _n);
}


void IGLCore::PolyharmonicDeformation()
{
	double z_max = 1.0;
	double z_dir = -0.03;
	int k = 2;
	bool resolve = true;
	Eigen::MatrixXd V, U;
	Eigen::VectorXd Z;
	Eigen::MatrixXi F;
	Eigen::VectorXi b;
	Eigen::VectorXd bc;

	using namespace Eigen;
	using namespace std;
	igl::readOBJ(p_readfilepath + "bump-domain.obj", V, F);
	U = V;
	// Find boundary vertices outside annulus
	typedef Matrix<bool, Dynamic, 1> VectorXb;
	VectorXb is_outer = (V.rowwise().norm().array() - 1.0) > -1e-15;
	VectorXb is_inner = (V.rowwise().norm().array() - 0.15) < 1e-15;
	VectorXb in_b = is_outer.array() || is_inner.array();
	igl::colon<int>(0, V.rows() - 1, b);
	b.conservativeResize(stable_partition(b.data(), b.data() + b.size(),
		[&in_b](int i)->bool {return in_b(i);}) - b.data());
	bc.resize(b.size(), 1);
	for (int bi = 0;bi < b.size();bi++)
	{
		bc(bi) = (is_outer(b(bi)) ? 0.0 : 1.0);
	}

	igl::harmonic(V, F, b, bc, k, Z);
	U.col(2) = z_max*Z;
	_v = U;
	_f = F;

	igl::per_face_normals(_v, _f, _n);
	igl::writeSTL(p_writefilepath + "Polyharmonic.stl", _v, _f, _n, true);


}

void IGLCore::BiharmonicDeformation(std::string dmatfile)
{
	Eigen::MatrixXd U, U_bc;
	Eigen::VectorXd Z;
	Eigen::VectorXi b;

	using namespace Eigen;
	using namespace std;
	U = _v;
	VectorXi S;
	igl::readDMAT(p_readfilepath + dmatfile, S);

	igl::colon<int>(0, _v.rows() - 1, b);
	b.conservativeResize(stable_partition(b.data(), b.data() + b.size(), [&S](int i)->bool {return S(i) >= 0;}) - b.data());

	// Boundary conditions directly on deformed positions
	U_bc.resize(b.size(), _v.cols());
	for (int bi = 0;bi < b.size();bi++)
	{
		//Translate
		U_bc.row(bi).x() = _v.row(b(bi)).x();
		U_bc.row(bi).y() = _v.row(b(bi)).y();
		U_bc.row(bi).z() = _v.row(b(bi)).z();
	}
	igl::harmonic(_v, _f, b, U_bc, 2., U);
	_v = U;
}

void IGLCore::MarchingCube()
{
	using namespace Eigen;
	using namespace igl;
	// number of vertices on the largest side
	const int s = 50;
	const RowVector3d Vmin = _v.colwise().minCoeff();
	const RowVector3d Vmax = _v.colwise().maxCoeff();
	const double h = (Vmax - Vmin).maxCoeff() / (double)s;
	const RowVector3i res = (s*((Vmax - Vmin) / (Vmax - Vmin).maxCoeff())).cast<int>();
	// create grid
	std::cout << "Creating grid..." << std::endl;
 	MatrixXd GV(res(0)*res(1)*res(2), 3);
	for (int zi = 0;zi < res(2);zi++)
	{
		const auto lerp = [&](const int di, const int d)->double
		{return Vmin(d) + (double)di / (double)(res(d) - 1)*(Vmax(d) - Vmin(d));};
		const double z = lerp(zi, 2);
		for (int yi = 0;yi < res(1);yi++)
		{
			const double y = lerp(yi, 1);
			for (int xi = 0;xi < res(0);xi++)
			{
				const double x = lerp(xi, 0);
				GV.row(xi + res(0)*(yi + res(1)*zi)) = RowVector3d(x, y, z);
			}
		}
	}
	// compute values
	std::cout << "Computing distances..." << std::endl;
	VectorXd S, B;
	VectorXi I;
	MatrixXd C, N;
	signed_distance(GV, _v, _f, SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C, N);
	// Convert distances to binary inside-outside data --> aliasing artifacts
	B = S;
	std::for_each(B.data(), B.data() + B.size(), [](double& b) {b = (b > 0 ? 1 : (b < 0 ? -1 : 0));});

	std::cout << "Marching cubes..." << std::endl;
	MatrixXd SV, BV;
	MatrixXi SF, BF;
	MatrixXd SN, BN;
	igl::copyleft::marching_cubes(S, GV, res(0), res(1), res(2), SV, SF);
	igl::copyleft::marching_cubes(B, GV, res(0), res(1), res(2), BV, BF);
	std::cout << "Marching cubes Complete" << std::endl;

	igl::per_face_normals(SV, SF, SN);
	igl::writeSTL(p_writefilepath + "Marching_Cube_result1.stl", SV, SF, SN);
	igl::per_face_normals(BV, BF, BN);
	igl::writeSTL(p_writefilepath + "Marching_Cube_result2.stl", BV, BF, BN);
	//igl::writeSTL(p_writefilepath + "Marching_Cube_result_origin.stl", V, F, Nor);

	//
	//using namespace Eigen;
	//using namespace std;
	//using namespace igl;
	//MatrixXi F, FF;
	//MatrixXd V, VV;
	//MatrixXi Nor, NorNor;
	//// Read in inputs as double precision floating point meshes
	//igl::readSTL(p_readfilepath + "LowerJawScan.stl", VV, FF, NorNor);
	//Eigen::MatrixXd SI;
	//Eigen::MatrixXd SVJ;
	//igl::remove_duplicate_vertices(VV, 0.000001, V, SI, SVJ);
	//F.resizeLike(FF);
	//for (int f = 0;f < FF.rows();f++)
	//{
	//	for (int c = 0;c < FF.cols();c++)
	//	{
	//		F(f, c) = SVJ(FF(f, c));
	//	}
	//}
	////igl::per_face_normals_stable(V, F, Nor);
	//// number of vertices on the largest side
	//const RowVector3d Vmin = V.colwise().minCoeff();
	//const RowVector3d Vmax = V.colwise().maxCoeff();
	//// create grid
	//std::cout << "Creating grid..." << std::endl;
	//int grid_x_size = 50;
	//int grid_y_size = 50;
	//int grid_z_size = 50;
	//MatrixXd GV(grid_x_size * grid_y_size * grid_z_size, 3);
	//for (int zi = 0;zi < grid_z_size;zi++)
	//{
	//	const auto lerp = [&](const int di, const int d, const int Vd)->double
	//	{return Vmin(Vd) + (double)di / (double)(d - 1)*(Vmax(Vd) - Vmin(Vd));};
	//	const double z = lerp(zi, grid_z_size, 2);
	//	for (int yi = 0;yi < grid_y_size;yi++)
	//	{
	//		const double y = lerp(yi, grid_y_size, 1);
	//		for (int xi = 0;xi < grid_x_size;xi++)
	//		{
	//			const double x = lerp(xi, grid_x_size, 0);
	//			GV.row(xi + grid_x_size*(yi + grid_y_size*zi)) = RowVector3d(x, y, z);
	//		}
	//	}
	//}
	//// compute values
	//std::cout << "Computing distances..." << std::endl;
	//VectorXd S, B;
	//{
	//	VectorXi I;
	//	MatrixXd C, N;
	//	signed_distance(GV, V, F, SIGNED_DISTANCE_TYPE_PSEUDONORMAL, S, I, C, N);
	//	// Convert distances to binary inside-outside data --> aliasing artifacts
	//	B = S;
	//	for_each(B.data(), B.data() + B.size(), [](double& b) {b = (b > 0 ? 1 : (b < 0 ? -1 : 0));});
	//}
	//
	///*std::ofstream writeFile;
	//writeFile.open(p_writefilepath + "SDF_file_test.txt");
	//int boundary_count = 0;

	//if (writeFile.is_open())
	//{
	//	for (int i = 0; i < B.rows(); i++)
	//		writeFile << B(i) << std::endl;
	//}
	//writeFile.close();*/
	//std::cout << "Marching cubes..." << std::endl;
	//MatrixXd SV, BV;
	//MatrixXi SF, BF;
	//MatrixXd SN, BN;
	//igl::copyleft::marching_cubes(S, GV, grid_x_size, grid_y_size, grid_z_size, SV, SF);
	//igl::copyleft::marching_cubes(B, GV, grid_x_size, grid_y_size, grid_z_size, BV, BF);
	//std::cout << "Marching cubes Complete" << std::endl;

	//igl::per_face_normals(SV, SF, SN);
	//igl::writeSTL(p_writefilepath + "Marching_Cube_result1.stl", SV, SF, SN);
	//igl::per_face_normals(BV, BF, BN);
	//igl::writeSTL(p_writefilepath + "Marching_Cube_result2.stl", BV, BF, BN);
	//igl::writeSTL(p_writefilepath + "Marching_Cube_result_origin.stl", V, F, Nor);
	//igl::per_face_normals(SV, SF, _n);
	//igl::writeSTL(p_writefilepath + "Marching_Cube_result.stl", SV, SF, _n, true);

}


void IGLCore::ShapeUp()
{
	Eigen::MatrixXd VQC;
	Eigen::MatrixXi FQC;
	Eigen::MatrixXi NQC;
	Eigen::MatrixXi E;
	Eigen::MatrixXi FQCtri;
	Eigen::MatrixXd PQC0, PQC1, PQC2, PQC3;
	// Euclidean-regular quad mesh
	Eigen::MatrixXd VQCregular;
	Eigen::MatrixXi FQCtriregular;
	Eigen::MatrixXd PQC0regular, PQC1regular, PQC2regular, PQC3regular;

	igl::ShapeupData su_data;



	// Scale for visualizing the fields
	double global_scale; //TODO: not used

	using namespace Eigen;
	using namespace std;

	// Load a quad mesh
	igl::readOBJ("D:/hbs/Code/LibIGL_ver2/LibIGL_ver2/model/basic_cube_Sub_3.quad.obj", VQC, FQC);

	FQCtri.resize(2 * FQC.rows(), 3);
	FQCtri << FQC.col(0), FQC.col(1), FQC.col(2),
		FQC.col(2), FQC.col(3), FQC.col(0);
	igl::slice(VQC, FQC.col(0).eval(), 1, PQC0);
	igl::slice(VQC, FQC.col(1).eval(), 1, PQC1);
	igl::slice(VQC, FQC.col(2).eval(), 1, PQC2);
	igl::slice(VQC, FQC.col(3).eval(), 1, PQC3);

	E.resize(FQC.size(), 2);
	E.col(0) << FQC.col(0), FQC.col(1), FQC.col(2), FQC.col(3);
	E.col(1) << FQC.col(1), FQC.col(2), FQC.col(3), FQC.col(0);

	VectorXi b(1); b(0) = 0;  //setting the first vertex to be the same.

	VectorXd wShape = VectorXd::Constant(FQC.rows(), 1.0);
	VectorXd wSmooth = VectorXd::Constant(E.rows(), 1.0);
	MatrixXd bc(1, 3); bc << VQC.row(0);

	VectorXi array_of_fours = VectorXi::Constant(FQC.rows(), 4);
	igl::shapeup_projection_function localFunction(igl::shapeup_regular_face_projection);

	su_data.maxIterations = 1000; // 200
	shapeup_precomputation(VQC, array_of_fours, FQC, E, b, wShape, wSmooth, su_data);
	shapeup_solve(bc, localFunction, VQC, su_data, false, VQCregular);

	// Convert the planarized mesh to triangles
	igl::slice(VQCregular, FQC.col(0).eval(), 1, PQC0regular);
	igl::slice(VQCregular, FQC.col(1).eval(), 1, PQC1regular);
	igl::slice(VQCregular, FQC.col(2).eval(), 1, PQC2regular);
	igl::slice(VQCregular, FQC.col(3).eval(), 1, PQC3regular);

	// Assign a color to each quad that corresponds to its planarity
	VectorXd angleRegularity(FQC.rows());
	quadAngleRegularity(VQCregular, FQC, angleRegularity);

	Eigen::MatrixXd SI;
	Eigen::MatrixXd SVJ;
	igl::remove_duplicate_vertices(VQCregular, 0.000001, _v, SI, SVJ);
	_f.resizeLike(FQCtri);
	for (int f = 0;f < FQCtri.rows();f++)
	{
		for (int c = 0;c < FQCtri.cols();c++)
		{
			_f(f, c) = SVJ(FQCtri(f, c));
		}
	}
	igl::per_face_normals_stable(_v, _f, _n);

}
void IGLCore::quadAngleRegularity(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Q, Eigen::VectorXd& angleRegularity)
{
	angleRegularity.conservativeResize(Q.rows());
	angleRegularity.setZero();
	for (int i = 0;i < Q.rows();i++) {
		for (int j = 0;j < 4;j++) {
			Eigen::RowVectorXd v21 = (V.row(Q(i, j)) - V.row(Q(i, (j + 1) % 4))).normalized();
			Eigen::RowVectorXd v23 = (V.row(Q(i, (j + 2) % 4)) - V.row(Q(i, (j + 1) % 4))).normalized();

			angleRegularity(i) += (abs(acos(v21.dot(v23)) - igl::PI / 2.0) / (igl::PI / 2.0)) / 4.0;
		}
	}
}


//////////////////////////제대로 안됨 수정 보완 필요 /////////////////////////

void IGLCore::Tetrahedralization()
{
	//Eigen::MatrixXd V;
	//Eigen::MatrixXi F;
	//Eigen::MatrixXd B;

	//// Tetrahedralized interior
	//Eigen::MatrixXd TV;
	//Eigen::MatrixXi TT;
	//Eigen::MatrixXi TF;

	//using namespace Eigen;
	//using namespace std;

	//// Load a surface mesh
	//igl::readOFF("/fertility.off", V, F);

	//// Tetrahedralize the interior
	//igl::copyleft::tetgen::tetrahedralize(V, F, "coadsflB", TV, TT, TF);

	//// Compute barycenters
	//igl::barycenter(TV, TT, B);

	//double t = double((4) + 1) / 9.0;

	//VectorXd v = B.col(2).array() - B.col(2).minCoeff();
	//v /= v.col(0).maxCoeff();

	//vector<int> s;

	//for (unsigned i = 0; i < v.size();++i)
	//	if (v(i) < t)
	//		s.push_back(i);

	//MatrixXd V_temp(s.size() * 4, 3);
	//MatrixXi F_temp(s.size() * 4, 3);

	//for (unsigned i = 0; i < s.size();++i)
	//{
	//	V_temp.row(i * 4 + 0).x() = TV.row(TT(s[i], 0)).x();
	//	V_temp.row(i * 4 + 0).y() = TV.row(TT(s[i], 0)).y();
	//	V_temp.row(i * 4 + 0).z() = TV.row(TT(s[i], 0)).z();


	//	V_temp.row(i * 4 + 1).x() = TV.row(TT(s[i], 1)).x();
	//	V_temp.row(i * 4 + 1).y() = TV.row(TT(s[i], 1)).y();
	//	V_temp.row(i * 4 + 1).z() = TV.row(TT(s[i], 1)).z();

	//	V_temp.row(i * 4 + 2).x() = TV.row(TT(s[i], 2)).x();
	//	V_temp.row(i * 4 + 2).y() = TV.row(TT(s[i], 2)).y();
	//	V_temp.row(i * 4 + 2).z() = TV.row(TT(s[i], 2)).z();

	//	V_temp.row(i * 4 + 3).x() = TV.row(TT(s[i], 3)).x();
	//	V_temp.row(i * 4 + 3).y() = TV.row(TT(s[i], 3)).y();
	//	V_temp.row(i * 4 + 3).z() = TV.row(TT(s[i], 3)).z();
	//	/*V_temp.row(i * 4 + 0) = TV.row(TT(s[i], 0));
	//	V_temp.row(i * 4 + 1) = TV.row(TT(s[i], 1));
	//	V_temp.row(i * 4 + 2) = TV.row(TT(s[i], 2));
	//	V_temp.row(i * 4 + 3) = TV.row(TT(s[i], 3));*/
	//	F_temp.row(i * 4 + 0) << (i * 4) + 0, (i * 4) + 1, (i * 4) + 3;
	//	F_temp.row(i * 4 + 1) << (i * 4) + 0, (i * 4) + 2, (i * 4) + 1;
	//	F_temp.row(i * 4 + 2) << (i * 4) + 3, (i * 4) + 2, (i * 4) + 0;
	//	F_temp.row(i * 4 + 3) << (i * 4) + 1, (i * 4) + 2, (i * 4) + 3;
	//}

}
void IGLCore::BooleanMesh()
{

	//	Eigen::MatrixXd VA, VB, VC;
	//	Eigen::VectorXi J, I;
	//	Eigen::MatrixXi FA, FB, FC;
	//	Eigen::MatrixXd NA, NB, NC;
	//	igl::MeshBooleanType boolean_type(
	//		igl::MESH_BOOLEAN_TYPE_UNION);
	//
	//	const char * MESH_BOOLEAN_TYPE_NAMES[] =
	//	{
	//		"Union",
	//		"Intersect",
	//		"Minus",
	//		"XOR",
	//		"Resolve",
	//	};
	//
	//
	//	using namespace Eigen;
	//	using namespace std;
	//	igl::readSTL("D:/hbs/Code/LibIGL_ver2/LibIGL_ver2/model/basic_cube.stl", VA, FA, NA);
	//	igl::readSTL("D:/hbs/Code/LibIGL_ver2/LibIGL_ver2/model/basic_cube_Sub_3.stl", VB, FB, NB);
	//	// Plot the mesh with pseudocolors
	//	
	//	// Initialize
	//	igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, boolean_type, VC, FC, J);
	//	
	//	igl::per_face_normals(VC, FC, NC);
	//	igl::writeSTL("D:/hbs/Code/LibIGL_ver2/LibIGL_ver2/model/result/boolean_test.stl", VC, FC, NC, true);
	//	
	//	cout <<
	//		"Press '.' to switch to next boolean operation type." << endl <<
	//		"Press ',' to switch to previous boolean operation type." << endl <<
	//		"Press ']' to push near cutting plane away from camera." << endl <<
	//		"Press '[' to pull near cutting plane closer to camera." << endl <<
	//		"Hint: investigate _inside_ the model to see orientation changes." << endl;
	//	
}

void IGLCore::Slice()
{
	using namespace Eigen;
	using namespace std;

	Eigen::MatrixXd V, BC;
	Eigen::VectorXd W;
	Eigen::MatrixXi T, F, G;
	double slice_z = 0.97;
	enum OverLayType
	{
		OVERLAY_NONE = 0,
		OVERLAY_INPUT = 1,
		OVERLAY_OUTPUT = 2,
		NUM_OVERLAY = 3,
	} overlay = OVERLAY_NONE;

	cout << "Usage:" << endl;
	cout << "[space]  toggle showing input mesh, output mesh or slice " << endl;
	cout << "         through tet-mesh of convex hull." << endl;
	cout << "'.'/','  push back/pull forward slicing plane." << endl;
	cout << endl;

	// Load mesh: (V,T) tet-mesh of convex hull, F contains facets of input
	// surface mesh _after_ self-intersection resolution
	igl::readMESH("model/big-sigcat.mesh", V, T, F);

	// Compute barycenters of all tets
	igl::barycenter(V, T, BC);

	// Compute generalized winding number at all barycenters
	cout << "Computing winding number over all " << T.rows() << " tets..." << endl;
	cout << V.rows() << " " << V.cols() << endl;
	cout << T.rows() << " " << T.cols() << endl;
	cout << F.rows() << " " << F.cols() << endl;
	cout << BC.rows() << " " << BC.cols() << endl;
	//cout << _v.rows() << _v.cols() << endl;
	cout << W.rows() << " " << W.cols() << endl;
	igl::winding_number(V, F, BC, W);

	cout << W.rows() << " " << W.cols() << endl;
	// Extract interior tets
	MatrixXi CT((W.array() > 0.5).count(), 4);
	{
		size_t k = 0;
		for (size_t t = 0;t < T.rows();t++)
		{
			if (W(t) > 0.5)
			{
				CT.row(k).x() = T.row(t).x();
				CT.row(k).y() = T.row(t).y();
				CT.row(k).z() = T.row(t).z();
				k++;
			}
		}
	}
	// find bounary facets of interior tets
	igl::boundary_facets(CT, G);
	// boundary_facets seems to be reversed...
	G = G.rowwise().reverse().eval();

	// normalize
	W = (W.array() - W.minCoeff()) / (W.maxCoeff() - W.minCoeff());

	// Plot the generated mesh
	//igl::opengl::glfw::Viewer viewer;

	//using namespace Eigen;
	//using namespace std;
	Eigen::Vector4d plane(
		0, 0, 1, -((1 - slice_z)*V.col(2).minCoeff() + slice_z*V.col(2).maxCoeff()));
	MatrixXd V_vis;
	MatrixXi F_vis;
	VectorXi J;
	{
		SparseMatrix<double> bary;
		// Value of plane's implicit function at all vertices
		const VectorXd IV =
			(V.col(0)*plane(0) +
				V.col(1)*plane(1) +
				V.col(2)*plane(2)).array()
			+ plane(3);
		igl::marching_tets(V, T, IV, V_vis, F_vis, J, bary);
	}
	VectorXd W_vis;
	cout << W.rows() << " " << W.cols() << endl;
	cout << J.rows() << " " << J.cols() << endl;
	cout << "J(0) = " << J[0] << endl;
	igl::slice(W, J, W_vis);
	cout << W_vis.rows() << " " << W_vis.cols() << endl;
	MatrixXd C_vis;
	// color without normalizing
	//igl::parula(W_vis, false, C_vis);
	cout << "Middle Point complete..." << endl;


	const auto & append_mesh = [&C_vis, &F_vis, &V_vis](
		const Eigen::MatrixXd & V,
		const Eigen::MatrixXi & F,
		const RowVector3d & color)
	{
		F_vis.conservativeResize(F_vis.rows() + F.rows(), 3);
		F_vis.bottomRows(F.rows()) = F.array() + V_vis.rows();
		V_vis.conservativeResize(V_vis.rows() + V.rows(), 3);
		V_vis.bottomRows(V.rows()) = V;
		C_vis.conservativeResize(C_vis.rows() + F.rows(), 3);
		C_vis.bottomRows(F.rows()).rowwise() = color;
	};
	switch (overlay)
	{
	case OVERLAY_INPUT:
		append_mesh(V, F, RowVector3d(1., 0.894, 0.227));
		break;
	case OVERLAY_OUTPUT:
		append_mesh(V, G, RowVector3d(0.8, 0.8, 0.8));
		break;
	default:
		break;
	}

	_v = W_vis;
	_f = F;
	//_n.setConstant();
	igl::per_face_normals(_v, _f, _n);
	cout << "Slice Function complete..." << endl;
	WriteSTL("Slice_result.stl", true);

}


void IGLCore::TestBed()
{
	using namespace Eigen;
	using namespace igl;
	MatrixXi F, FF;
	MatrixXd V, VV;
	MatrixXi N, NN;
	// Read in inputs as double precision floating point meshes
	igl::readSTL(p_readfilepath + "basic_cube.stl", VV, FF, NN);
	Eigen::MatrixXd SI;
	Eigen::MatrixXd SVJ;
	igl::remove_duplicate_vertices(VV, 0.000001, V, SI, SVJ);
	F.resizeLike(FF);
	for (int f = 0;f < FF.rows();f++)
	{
		for (int c = 0;c < FF.cols();c++)
		{
			F(f, c) = SVJ(FF(f, c));
		}
	}
	igl::per_face_normals_stable(V, F, N);
	std::cout << V << std::endl;

	V.conservativeResize(3, 8);
	std::cout << "###############" << std::endl;
	std::cout << V << std::endl;
}