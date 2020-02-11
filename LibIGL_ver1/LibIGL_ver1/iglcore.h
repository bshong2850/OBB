#pragma once
#ifndef CORE_H
#define CORE_H
#include <igl/readSTL.h>
#include <igl/writeSTL.h>

//duplication
#include <igl/remove_duplicate_vertices.h>
#include <igl/remove_duplicates.h>

//Compute Normal
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>

//Smoothing
#include <igl/barycenter.h>
#include <igl/cotmatrix.h>
#include <igl/doublearea.h>
#include <igl/grad.h>
#include <igl/jet.h>
#include <igl/massmatrix.h>

//Decimation
#include <igl/circulation.h>
#include <igl/collapse_edge.h>
#include <igl/edge_flaps.h>
#include <igl/shortest_edge_and_midpoint.h>
#include <igl/read_triangle_mesh.h>
#include <igl/decimate.h>
#include <Eigen/Core>
#include <set>

//Subdivision
#include <igl/upsample.h>
#include <igl/loop.h>

//ShapeUp
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/jet.h>
#include <igl/shapeup.h>
#include <igl/quad_planarity.h>
#include <igl/readDMAT.h>
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/PI.h>
#include <vector>
#include <cstdlib>

//Marching Cube
#include <igl/copyleft/marching_cubes.h>
#include <igl/signed_distance.h>
#include <igl/read_triangle_mesh.h>

//Tetrahedrailzation
//#include <igl/copyleft/tetgen/tetrahedralize.h>
//#include <igl/readOFF.h>
//#include <igl/barycenter.h>


//Slice
#include <igl/readMESH.h>
#include <igl/floor.h>
#include <igl/slice.h>
#include <igl/slice_into.h>
#include <igl/barycenter.h>
#include <igl/boundary_facets.h>
#include <igl/marching_tets.h>
#include <igl/winding_number.h>

//Deformation
#include <igl/colon.h>
#include <igl/harmonic.h>

//Deformation Sample
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readDMAT.h>
#include <algorithm>

#include <igl/polygon_mesh_to_triangle_mesh.h>



#include <iostream>
#include <windows.h>

class IGLCore
{
private:


	Eigen::MatrixXd p_originV;
	Eigen::MatrixXi p_originF;
	Eigen::MatrixXd p_originN;

	std::string p_readfilepath = "model/";
	std::string p_writefilepath = "model/result/";;
protected:

public:
	IGLCore();
	~IGLCore();

	Eigen::MatrixXd _v;
	Eigen::MatrixXi _f;
	Eigen::MatrixXd _n;

	void Initialize();
	void PrintInformation();

	void ReadFile(std::string filename, std::string extension);
	void ReadSTL(std::string filename, bool _duplicate = true);
	void ReadOBJ(std::string filename);
	void WriteDMAT(std::string filename);

	std::ofstream writeFile;
	void WriteSTL(std::string filename, bool _ascii = false);

	enum NormalType
	{
		FACE_NORMAL,
		VERTEX_NORMAL,
		CORNER_NORMAL
	};
	void ComputeNormal(int mode);

	void Smoothing();
	void Decimation();

	void Subdivision();

	void PolyharmonicDeformation();
	void BiharmonicDeformation(std::string dmatfile);

	void MarchingCube();
	void TestBed();

	void Tetrahedralization();
	//Quad Mesh¸¸ °¡´É
	void ShapeUp();
	void quadAngleRegularity(const Eigen::MatrixXd& V, const Eigen::MatrixXi& Q, Eigen::VectorXd& angleRegularity);

	void BooleanMesh();
	void Slice();


};



#endif