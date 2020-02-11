#pragma once
#ifndef CORE_H
#define CORE_H
#include <igl/readSTL.h>
#include <igl/writeSTL.h>

//duplication
#include <igl/remove_duplicate_vertices.h>

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

//Deformation
#include <igl/colon.h>
#include <igl/harmonic.h>

//Deformation Sample
#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/readDMAT.h>
#include <algorithm>



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


	void ReadSTL(std::string filename, bool _duplicate = true);
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

	void BiharmonicDeformation(std::string dmatfile);
	void PolyharmonicDeformation();
};



#endif