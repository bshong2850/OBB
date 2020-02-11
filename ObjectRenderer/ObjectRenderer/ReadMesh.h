#pragma once
#include <fstream>
#include<string> 
#include"Vector.h"
#include<vector>
#include<GL\glut.h>
#include<math.h>
#include<algorithm>

#include "STL_Loader.h"
 
using namespace std;
class face;

class face
{
public:
	int					f_Vertex[3];
	vec3				m_Normal;
	vec3				m_centroid;
	int					ID;
	bool				Select;
	bool				SmoothRegion;


public:
	face() {};
	face(int i, int j, int k) {
		f_Vertex[0] = i;
		f_Vertex[1] = j;
		f_Vertex[2] = k; 
		SmoothRegion = false;
		Select = false;
	}

	~face() {}
	std::vector<face*>		Neighbor_face;
};

class vertex
{
public:
	vertex() {}
	vertex(double x, double y, double z)
	{
		setVertex(x, y, z);
	}
	~vertex() {}

public:
	
	double			m_InvMass;


	vector<face*>		Neighbor_face;
	vector<vertex*> Neighbor_Vertex;
	vec3			pos;
	vec3			ori_pos;
	vec3			m_Vel;
	vec3			m_NewPos; 
	int id;
	bool Smooth;
	vec3			m_normal;
	vec2		textcoord;
	vec3		Color;

 

	void setVertex(double x, double y, double z)
	{
		pos.x = x;
		pos.y = y;
		pos.z = z;
		ori_pos.x = x;
		ori_pos.y = y;
		ori_pos.z = z;
		Smooth = false;
		Color = vec3(0.95, 0.88, 0.8);
	}

};

class edge
{
public: 
	~edge() {};
	int					e_Vertex[2];
	vector<face*>		Neighborface;
	int					e_id;
	double				e_Length;
	bool Smooth;

public:
	edge() { e_id = -1; }
	int p1, p2;
	edge(int _p1, int _p2)
	{
		p1 = _p1;
		p2 = _p2;
		if (p1 > p2)
		{
			int t = p1;
			p1 = p2;
			p2 = t;
		}
		e_Vertex[0] = p1;
		e_Vertex[1] = p2;
		Smooth = false;
	}
 
	bool operator <(edge& e)
	{
		if (p1 - e.p1) return p1 < e.p1;
		if (p2 - e.p2) return p2 < e.p2;
		return false;
	}
	bool operator ==(edge& e)
	{
		return (p1 == e.p1 && p2 == e.p2) || (p1 == e.p2 && p2 == e.p1);
	}
	bool operator !=(edge& e)
	{
		return (p1 != e.p1 && p2 != e.p2) || (p1 != e.p2 && p2 != e.p1);
	}
};

class mesh
{
public:
	vector<vertex>			m_Vertices;
	vector<face>			m_Faces;
	vector<edge>			m_Edges;

	std::string				stl_header;
	std::string				stl_name;
public:
	mesh(void) {};
	~mesh(void) {};
public:

	void	vertexId(void)
	{
		int number_v = (int)m_Vertices.size();
		for (int i = 0; i < number_v; i++)
		{
			vertex& v = m_Vertices[i];
			v.id = i;
		}
	}

	void	Adjacent_face(void)
	{
		int number_faces = (int)m_Faces.size();
		for (int i = 0; i < number_faces; i++)
		{
			face& f = m_Faces[i];
			for (int j = 0; j < 3; j++)
			{
				int ver_id = f.f_Vertex[j];
				vertex& v = m_Vertices[ver_id];
				v.Neighbor_face.push_back(&f);
			}
		}

		for (int i = 0; i < m_Vertices.size(); i++)
		{
			vertex &v = m_Vertices[i];
			vec3 normaltotal = vec3(0, 0, 0);

			for (int j = 0; j < v.Neighbor_face.size(); j++)
			{
				face* f = v.Neighbor_face[j];
				for (int k = 0; k < 3; k++)
				{
					normaltotal = normaltotal + f->m_Normal;
				}
			}
			v.m_normal = normaltotal / v.Neighbor_face.size();
			v.m_normal.Normalize();
		}

		for (int i = 0; i < m_Vertices.size(); i++)
		{
			vertex& v = m_Vertices[i];

			for (int j = 0; j < v.Neighbor_face.size(); j++)
			{
				face* f = v.Neighbor_face[j];

				for (int k = 0; k < 3; k++)
				{
					if (f->f_Vertex[k] != v.id)
						v.Neighbor_Vertex.push_back(&m_Vertices[f->f_Vertex[k]]);
				}
			}
			sort(v.Neighbor_Vertex.begin(), v.Neighbor_Vertex.end());
			v.Neighbor_Vertex.erase(unique(v.Neighbor_Vertex.begin(), v.Neighbor_Vertex.end()), v.Neighbor_Vertex.end());

		}

	}

	void	edge_Adjacent_face(void)
	{
		int number_edges = (int)m_Edges.size();
		for (int i = 0; i < number_edges; i++)
		{
			edge& e = m_Edges[i];
			int		ver_id1 = e.e_Vertex[0];
			int		ver_id2 = e.e_Vertex[1];
			vertex& v1 = m_Vertices[ver_id1];
			vertex& v2 = m_Vertices[ver_id2];
			//e length
			e.e_Length = v1.pos.dist(v2.pos);
			int		num_neighbor1 = v1.Neighbor_face.size();
			int		num_neighbor2 = v2.Neighbor_face.size();

			for (int j = 0; j < num_neighbor2; j++)
			{
				for (int k = 0; k < num_neighbor1; k++)
				{
					face* f = v2.Neighbor_face[j];
					face* f1 = v1.Neighbor_face[k];
					if (f1 == f)
					{
						e.Neighborface.push_back(f);
						e.Neighborface.push_back(f1);
					}
				}
			}
			e.Neighborface.erase(unique(e.Neighborface.begin(), e.Neighborface.end()), e.Neighborface.end());
		}
	}

	void	Compute_Face_Normal(void)
	{
		int number_faces = (int)m_Faces.size();
		for (int i = 0; i < number_faces; i++)
		{
			face	&f = m_Faces[i];
			f.Select = false;
			vec3	v0 = m_Vertices[f.f_Vertex[0]].pos;
			vec3	v1 = m_Vertices[f.f_Vertex[1]].pos;
			vec3	v2 = m_Vertices[f.f_Vertex[2]].pos;
			vec3	va = v1 - v0;
			vec3	vb = v2 - v0;
			vec3	vc = va.Cross(vb);
			vec3	vd = va.Cross(vb);
			vc.Normalize();
			f.m_Normal = vc;								// compute normal
			f.m_centroid = (v0 + v1 + v2) / 3;				//compute centroid 
			f.ID = i;
		}
	}
	void	Compute_partFace_Normal(void)
	{
		int number_faces = (int)m_Faces.size();
		for (int i = 0; i < number_faces; i++)
		{

			face& f = m_Faces[i];
			if (f.Select == true)
			{
				vec3	v0 = m_Vertices[f.f_Vertex[0]].pos;
				vec3	v1 = m_Vertices[f.f_Vertex[1]].pos;
				vec3	v2 = m_Vertices[f.f_Vertex[2]].pos;
				vec3	va = v1 - v0;
				vec3	vb = v2 - v0;
				vec3	vc = va.Cross(vb);
				vec3	vd = va.Cross(vb);
				vc.Normalize();
				f.m_Normal = vc;								// compute normal
				f.m_centroid = (v0 + v1 + v2) / 3;				//compute centroid 
																//	f.ID = i;
			}

		}
	}


	vec3   Compute_Min()
	{
		vec3			min;
		double			minY;
		double			minV = DBL_MAX;
		double			minVerticesY;
		double			minX;
		double			min_X_V = DBL_MAX;
		double			minVerticesX;
		double			minz;
		double			min_z_V = DBL_MAX;
		double			minVerticesz;
		int				number_f = (int)m_Faces.size();
		for (int i = 0; i < number_f; i++)
		{
			face& f = m_Faces[i];


			double verticeY = f.m_centroid.y;
			if (verticeY < minV)
			{
				minV = verticeY;
				minVerticesY = verticeY;
			}
			double verticeX = f.m_centroid.x;
			if (verticeX < min_X_V)
			{
				min_X_V = verticeX;
				minVerticesX = verticeX;
			}
			double verticez = f.m_centroid.z;
			if (verticez < min_z_V)
			{
				min_z_V = verticez;
				minVerticesz = verticez;
			}

		}
		minY = minVerticesY;
		minX = minVerticesX;
		minz = minVerticesz;
		min = vec3(minX, minY, minz);
		return min;
	}



	vec3   Compute_Max()
	{
		vec3			max;
		double			maxY;
		double			maxV = DBL_MIN;
		double			maxVerticesY;
		double			maxX;
		double			max_X_V = DBL_MIN;
		double			maxVerticesX;
		double			maxz;
		double			max_z_V = DBL_MIN;
		double			maxVerticesz;
		int				number_f = (int)m_Faces.size();
		for (int i = 0; i < number_f; i++)
		{
			face& f = m_Faces[i];

			double verticeY = f.m_centroid.y;
			if (verticeY > maxV)
			{
				maxV = verticeY;
				maxVerticesY = verticeY;
			}
			double verticeX = f.m_centroid.x;
			if (verticeX > max_X_V)
			{
				max_X_V = verticeX;
				maxVerticesX = verticeX;
			}
			double verticez = f.m_centroid.z;
			if (verticez > max_z_V)
			{
				max_z_V = verticez;
				maxVerticesz = verticez;
			}
		}

		maxY = maxVerticesY;
		maxX = maxVerticesX;
		maxz = maxVerticesz;
		max = vec3(maxX, maxY, maxz);
		return max;
	}

	void TextcoordinatesInit()
	{
		vec3 min = Compute_Min();
		vec3 max = Compute_Max();
		float x_size = max.x - min.x;
		float y_size = max.y - min.y;
		float z_size = max.z - min.z;
		cout << min.x << endl;

		for (int i = 0; i < m_Vertices.size(); i++)
		{
			vertex& v = m_Vertices[i];
			//	v.textcoord = vec2(v.pos.x*182 / m_Vertices.size(), v.pos.y*200 / m_Vertices.size());
			// v.textcoord = vec2(v.pos.x*0.8 / x_size, v.pos.y*0.65/y_size) ;
			v.textcoord = vec2(v.pos.x*0.8 / x_size, v.pos.y * 0.67 / y_size);

		}
	}

	void	Render_Solid(void)
	{
		int number_faces = (int)m_Faces.size();
		for (int i = 0; i < number_faces; i++)
		{
			face &f = m_Faces[i];
			if (f.SmoothRegion == false)
				glColor3f(0.95, 0.88, 0.8);
			else if (f.SmoothRegion == true)
				glColor3f(1.0, 0.0, 0.0);

			glBegin(GL_POLYGON);
			for (int j = 0; j < 3; j++)
			{
				glNormal3d(f.m_Normal.x, f.m_Normal.y, f.m_Normal.z);
				glVertex3f(m_Vertices[f.f_Vertex[j]].pos.x, m_Vertices[f.f_Vertex[j]].pos.y, m_Vertices[f.f_Vertex[j]].pos.z);
			}
			glEnd();
		}

	}

	void SelectDesicion(void)
	{
		int number_faces = (int)m_Faces.size();
		for (int i = 0; i < number_faces; i++)
		{
			face& f = m_Faces[i];
			if (f.SmoothRegion == true)
			{
				vertex& v1 = m_Vertices[f.f_Vertex[0]];
				vertex& v2 = m_Vertices[f.f_Vertex[1]];
				vertex& v3 = m_Vertices[f.f_Vertex[2]];
				v1.Smooth = true;
				v2.Smooth = true;
				v3.Smooth = true;

			}
		}

		int number_edges = (int)m_Edges.size();
		for (int i = 0; i < number_edges; i++)
		{
			edge& e = m_Edges[i];
			int		ver_id1 = e.e_Vertex[0];
			int		ver_id2 = e.e_Vertex[1];
			vertex& v1 = m_Vertices[ver_id1];
			vertex& v2 = m_Vertices[ver_id2];
			if (v1.Smooth == true || v2.Smooth == true)
				e.Smooth = true;
		}

	}
	vec3 rotate(vec3 a, vec3 ori, vec3 theta)
	{
		a = a - ori;
		double A = cos(theta.x);
		double B = sin(theta.x);
		double C = cos(theta.y);
		double D = sin(theta.y);
		double E = cos(theta.z);
		double F = sin(theta.z);

		vec3 p1, p2, p3;
		float rot[4][4] = { C * E, -C * F, -D, 0,
			-B * D * E + A * F, B * D * F + A * E, -B * C, 0,
			A * D * E + B * F, -A * D * F + B * E, A * C, 0,
			0, 0, 0, 1
		};

		float b[4] = { a.x , a.y, a.z,0 };

		float c[4] = { 0,0,0,0 };
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				c[i] += rot[i][j] * b[j];
			}
		}
		a.x = c[0];
		a.y = c[1];
		a.z = c[2];

		//translate to the original position of the point
		a = a + ori;
		return a;
	}

	void RotationByAngle(vec3 angle)
	{
		for (int i = 0; i < m_Vertices.size(); i++)
		{
			vertex& v = m_Vertices[i];
			vec3 ori = vec3(23, 15, 15);
			v.pos = rotate(v.pos, ori, angle);
		}
	}

};


 
class ReadMesh
{
public:
	ReadMesh();
	~ReadMesh();
public:
	mesh				TriMesh;
public:
	void	ReadMeshOBJ(std::string filename);
	void	ReadMeshSTL(std::string filename);
	void	RenderSTL();

	std::ifstream		bvfile;
	void	ReadBV(std::string filename);
	void	RenderBV();
};

