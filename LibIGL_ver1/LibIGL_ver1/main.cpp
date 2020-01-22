#include "iglcore.h"
int main(int argc, char *argv[])
{
	IGLCore IC;
	IC.ReadSTL("basic_cube.stl", true);

	/*std::cout << IC._v << std::endl;
	std::cout << "#################################" << std::endl;*/
	//IC.ComputeNormal(IC.VERTEX_NORMAL);
	/*for (int i = 0; i<3; i++)
		IC.Smoothing();*/
	for (int i = 1; i <= 100000; i++)
	{
		IC.BiharmonicDeformation();
		if (i % 1000 == 0)
			IC.WriteSTL("BiharmonicDeformation_basic_cube_" + std::to_string(i) + ".stl", true);

	}
	/*std::cout << IC._v << std::endl;
	std::cout << "#################################" << std::endl;*/
	return 0;

}