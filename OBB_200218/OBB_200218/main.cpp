#include "OBB.h"
#include <iostream>

int main() {

	SHAPE a, b;
	a.Setparameter();
	b.Setparameter();

	OBB obb;
	if (obb.ComputeOBB(a, b))
		std::cout << "crash" << std::endl;
	else
		std::cout << "Not crash" << std::endl;

}