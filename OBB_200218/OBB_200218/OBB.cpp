
#include "OBB.h"
#include <iostream>
using namespace std;

SHAPE::SHAPE()
{
	p_centery = 0;
	p_centerx = 0;
	p_height = 0;
	p_width = 0;
	p_rot = 0;
}

SHAPE::SHAPE(double centery, double centerx, double height, double width, double rot)
{
	p_centery = centery;
	p_centerx = centerx;
	p_height = height;
	p_width = width;
	p_rot = rot;
}

SHAPE::~SHAPE()
{

}

void SHAPE::Setparameter()
{
	double centery = 0;
	double centerx = 0;
	double height = 0;
	double width = 0;
	double rot = 0;
	std::cout << "Data Setting" << std::endl;
	std::cout << "center X Coordinate = ";
	std::cin >> centerx;
	std::cout << "center y Coordinate = ";
	std::cin >> centery;
	std::cout << "height distance = ";
	std::cin >> height;
	std::cout << "width distance = ";
	std::cin >> width;
	std::cout << "rotation value = ";
	std::cin >> rot;
}


VECTOR SHAPE::addVector(VECTOR a, VECTOR b) { //vector plus
	VECTOR ret;
	ret.x = a.x + b.x;
	ret.y = a.y + b.y;
	return ret;
}

double SHAPE::absDotVector(VECTOR a, VECTOR b) { //vector inner
	return abs(a.x * b.x + a.y * b.y);
}

double SHAPE::Deg2Rad(double deg) { //deg -> rad
	return deg / 180 * 3.141592;
}

VECTOR SHAPE::getDistanceVector(SHAPE a, SHAPE b) { //distance vector
	VECTOR ret;
	ret.x = (a.p_centerx + a.p_width / 2) - (b.p_centerx + b.p_width / 2);
	ret.y = (a.p_centery + a.p_height / 2) - (b.p_centery + b.p_height / 2);
	return ret;
}

VECTOR SHAPE::getHeightVector(SHAPE a) { //height vector
	VECTOR ret;
	ret.x = a.p_height * cos(Deg2Rad(a.p_rot - 90)) / 2;
	ret.y = a.p_height * sin(Deg2Rad(a.p_rot - 90)) / 2;
	return ret;
}

VECTOR SHAPE::getWidthVector(SHAPE a) { //width vector
	VECTOR ret;
	ret.x = a.p_width * cos(Deg2Rad(a.p_rot)) / 2;
	ret.y = a.p_width * sin(Deg2Rad(a.p_rot)) / 2;
	return ret;
}

VECTOR SHAPE::getUnitVector(VECTOR a) { //unit vector
	VECTOR ret;
	double size;
	size = sqrt(pow(a.x, 2) + pow(a.y, 2));
	ret.x = a.x / size;
	ret.y = a.y / size;
	return ret;
}

bool OBB::ComputeOBB(SHAPE a, SHAPE b) { //final check
	VECTOR dist = getDistanceVector(a, b);
	VECTOR vec[4] = { getHeightVector(a), getHeightVector(b), getWidthVector(a), getWidthVector(b) };
	VECTOR unit;
	for (int i = 0; i < 4; i++) {
		double sum = 0;
		unit = getUnitVector(vec[i]);
		for (int j = 0; j < 4; j++) {
			sum += absDotVector(vec[j], unit);
		}
		if (absDotVector(dist, unit) > sum) {
			return false;
		}
	}
	return true;
}
