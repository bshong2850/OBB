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
	IC.ReadSTL("LowerJawScan.stl", true);
	
	IC.WriteDMAT("LowerJawScan.dmat");
	IC.BiharmonicDeformation("LowerJawScan.dmat");
	IC.WriteSTL("LowerJawScan_deform_first_zero_1.stl", true);
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