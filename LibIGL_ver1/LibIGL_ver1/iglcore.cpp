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

void IGLCore::ReadSTL(std::string filename, bool _duplicate)
{
	igl::readSTL(p_readfilepath + filename, p_originV, p_originF, p_originN);
	if (_duplicate)
	{
		Eigen::MatrixXi SVI;
		Eigen::MatrixXi SVJ;
		igl::remove_duplicate_vertices(p_originV, p_originF, 0.0000001, _v, SVI, SVJ, _f);
		_n = p_originN;
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

void IGLCore::BiharmonicDeformation()
{
	using namespace Eigen;
	using namespace std;
	int k = 1;

	Eigen::VectorXd Z;
	Eigen::VectorXi b;
	Eigen::SparseMatrix<double> L;
	Eigen::VectorXd bc;

	// Find boundary vertices outside annulus
	typedef Matrix<bool, Dynamic, 1> VectorXb;
	VectorXb is_outer = (_v.rowwise().norm().array() - 1.0) > -1e-15;
	VectorXb is_inner = (_v.rowwise().norm().array() - 0.15) < 1e-15;
	VectorXb in_b = is_outer.array() || is_inner.array();
	igl::colon<int>(0, _v.rows() - 1, b);
	b.conservativeResize(stable_partition(b.data(), b.data() + b.size(),
		[&in_b](int i)->bool {return in_b(i);}) - b.data());
	bc.resize(b.size(), 1);
	for (int bi = 0;bi < b.size();bi++)
	{
		bc(bi) = (is_outer(b(bi)) ? 0.0 : 1.0);
	}


	//igl::harmonic(_v, _f, b, bc, k, Z);
	igl::harmonic(_v, _f, k, L);
}