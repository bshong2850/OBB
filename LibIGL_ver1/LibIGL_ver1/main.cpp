#include "iglcore.h"
int main(int argc, char *argv[])
{
	IGLCore IC;
	/*C.ReadSTL("basic_cube.stl", true);

	std::cout << IC._v << std::endl;
	std::cout << "#################################" << std::endl;
	std::cout << IC._v.row(0) << std::endl << IC._v.row(1) << std::endl;

	IC._v.row(0).x() = IC._v.row(1).x();
	IC._v.row(0).y() = IC._v.row(1).y();
	IC._v.row(0).z() = IC._v.row(1).z();
	std::cout << IC._v.row(0) << std::endl << IC._v.row(1) << std::endl;*/
	//std::string filename = "LowerJawScan";
	//std::string filename = "LowerJawScan";
	//IC.BooleanMesh();
	//IC.WriteSTL("test_aslkdfj.stl", true);
	//IC.WriteSTL(filename + "_decimation.stl", true);
	//std::cout << i << std::endl;


	std::string filename = "LowerJawScan";
	std::string extension = "stl";
	IC.ReadFile(filename, extension);;
	/*std::string filename = "armadillo";
	std::string extension = "obj";
	IC.ReadFile(filename, extension);*/
	//IC.ShapeUp();
	//IC.WriteSTL(filename + "_ShapeUp_1.stl", false);
	IC.Subdivision();
	IC.WriteSTL("Subdivision_" + filename + ".stl", false);
	//IC.TestBed();
	//IC.PolyharmonicDeformation();

	//IC.Decimation();
	//IC.WriteSTL(filename + "_decimation.stl", false);
	//IC.Slice();
	//std::string filename = "LowerJawScan";
	/*IC.ReadSTL(filename + ".stl", true);
	
	IC.WriteDMAT(filename + ".dmat");
	IC.BiharmonicDeformation(filename + ".dmat");
	IC.WriteSTL(filename + "_deform_first_zero_1.stl", true);*/
	//IC.WriteSTL("BiharmonicDeformation_result.stl", true);
	//IC.ComputeNormal(IC.VERTEX_NORMAL);
	/*for (int i = 0; i<3; i++)
		IC.Smoothing();*/
		/*for (int i = 1; i <= 100000; i++)
		{
			IC.BiharmonicDeformation();
			if (i % 1000 == 0)
				IC.WriteSTL("BiharmonicDeformation_basic_cube_" + std::to_string(i) + ".stl", true);

		}*/
	/*std::cout << IC._v << std::endl;
	std::cout << "#################################" << std::endl;*/
	return 0;

}