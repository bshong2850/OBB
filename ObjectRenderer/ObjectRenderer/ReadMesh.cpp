#define STLLOADER_IMPLEMENTATION
#include "ReadMesh.h"


ReadMesh::ReadMesh()
{
}

ReadMesh::~ReadMesh()
{
	bvfile.close();
}

void ReadMesh::ReadMeshOBJ(std::string filename)
{
	TriMesh.m_Faces.clear();
	TriMesh.m_Edges.clear();
	TriMesh.m_Vertices.clear();
	ifstream file;
	file.open(filename);
	double x, y, z;
	int i, j, k;
	if (file.fail())
	{
		cerr << "Cannot open file: " << filename << " exiting." << endl;
		return;
	}
	while (!file.eof())
	{
		char n;
		file >> n;
		string line;
		if (n == '#')
		{
			getline(file, line);
		}
		else if (n == 'v')
		{
			file >> x >> y >> z;
			vertex v(x, y, z);
			v.m_InvMass = 0.1;
			v.m_normal = vec3(0, 0, 0);
			v.m_Vel = vec3(0.0, 0.0, 0.0);

			TriMesh.m_Vertices.push_back(v);
		}
		else if (n == 'f')
		{
			file >> i >> j >> k;
			i--; j--; k--;

			TriMesh.m_Faces.push_back(face(i, j, k));
			TriMesh.m_Edges.push_back(edge(i, j));
			TriMesh.m_Edges.push_back(edge(k, j));
			TriMesh.m_Edges.push_back(edge(i, k));
		}
	}

	sort(TriMesh.m_Edges.begin(), TriMesh.m_Edges.end());
	TriMesh.m_Edges.erase(unique(TriMesh.m_Edges.begin(), TriMesh.m_Edges.end()), TriMesh.m_Edges.end());
	file.close();

	TriMesh.vertexId();
	TriMesh.Compute_Face_Normal();
	TriMesh.Adjacent_face();
	TriMesh.edge_Adjacent_face();
	TriMesh.TextcoordinatesInit();
}

stlloader::Mesh stlmesh;
void ReadMesh::ReadMeshSTL(std::string filename)
{
	stlloader::parse_file(filename, stlmesh);
	return ;
}

void ReadMesh::RenderSTL()
{
	for (auto & facet : stlmesh.facets)
	{
		glColor3f(0.95, 0.88, 0.8);

		glBegin(GL_POLYGON);
		for (int j = 0; j < 3; j++)
		{
			glNormal3d(facet.normal.x, facet.normal.y, facet.normal.z);
			glVertex3f(facet.vertices[j].x, facet.vertices[j].y, facet.vertices[j].z);
		}
		glEnd();

	}


	//stlloader::print(mesh);

	return;
}


void ReadMesh::ReadBV(std::string filename)
{
	bvfile;
	bvfile.open(filename);
	if (bvfile.fail())
	{
		cerr << "Cannot open file" << endl;
		return;
	}
	double x, y, z;
	while (!bvfile.eof())
	{
		bvfile >> x >> y >> z;
		vertex v(x, y, z);
		TriMesh.m_Vertices.push_back(v);
	}
	bvfile.close();
	return;
}

void ReadMesh::RenderBV()
{
	for (int i=0; i<TriMesh.m_Vertices.size(); i++)
	{
		//glColor3f(0.95, 0.88, 0.8);
		glColor3f(1,1,1);
		glPointSize(1.0f);
		glBegin(GL_POINTS);
		glNormal3d(0, 1, 0);
		glVertex3f(TriMesh.m_Vertices[i].pos.x,
			TriMesh.m_Vertices[i].pos.y,
			TriMesh.m_Vertices[i].pos.z);
		glEnd();
	}
	//stlloader::print(mesh);

	return;
}

