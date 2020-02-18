#pragma once

#ifndef OBB_H
#define OBB_H

struct VECTOR { //define vector
	double x, y;
};

class SHAPE
{
private:
	double p_centerx, p_centery, p_height, p_width, p_rot;

protected:

public:

	SHAPE();
	SHAPE(double top, double left, double height, double width, double rot);
	~SHAPE();

	void Setparameter();
	VECTOR addVector(VECTOR a, VECTOR b);
	double absDotVector(VECTOR a, VECTOR b);
	double Deg2Rad(double deg);
	VECTOR getDistanceVector(SHAPE a, SHAPE b);
	VECTOR getHeightVector(SHAPE a);
	VECTOR getWidthVector(SHAPE a);
	VECTOR getUnitVector(VECTOR a);
};

class OBB : public SHAPE
{
public:
	bool ComputeOBB(SHAPE a, SHAPE b);
};


#endif