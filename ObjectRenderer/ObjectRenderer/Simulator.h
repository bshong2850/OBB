

#pragma once

#include <iostream>
#include <vector>
#include <GL/glut.h>
#include <GL/gl.h>
#include <GL/glu.h>

#include "ReadMesh.h"
 

using namespace std;


class Simulator
{
public:
	Simulator(void);
	~Simulator(void);

public:
	void					Initialize(std::string filepath);
 	void					Update();
	void					Render();
	void					Lighting(void);
	void					DrawGround(void);

public:
	ReadMesh			readobj;
	std::string			extension;
	vec3				ground;
	float timsStep;
};

