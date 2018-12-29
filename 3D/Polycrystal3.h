#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"
#include "PolycrMesh.h"
#include "CrysesShell.h"

using std::ifstream;
using std::ofstream;
using std::unique_ptr;
using std::vector;
using std::string;

class Polycrystal3
{
private:
	double _preferredLength;

	PolycrMesh* lastTriangulation = nullptr;

	vector<Crystallite3*> crystallites;

	vector<ShellFacet3*> _shellFacets;
	vector<ShellEdge3*> _shellEdges;
	vector<ShellVertex3*> _shellVertexes;

	vector<unique_ptr<Facet3>*> _startFrontFacets;
	vector<unique_ptr<Edge3>*> _startFrontEdges;
	vector<unique_ptr<Vertex3>*> _startFrontVertexes;

	ShellEdge3* findShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const;
	unique_ptr<Edge3>* findStartFrontEdge(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1) const;

	template <class T>
	void removePtrsToNullptr(vector<unique_ptr<T>*>& vec);

	void setLinksWithShell();
	void removePtrsToNullptrFromVectors();
	void triangulateShell();
	void startFrontDelaunayPostprocessing();

	void inputData();

public:
	void generateMeshNoStruct(const double preferredLength);
	void generateMeshNoStruct(string filename, const double preferredLength);
	void generateMeshNoStruct(const CrysesShell& crysesShell, const double preferredLength);
	PolycrMesh* structurizeMesh();
	PolycrMesh* generateMesh(const double preferredLength);
	PolycrMesh* generateMesh(string filename, const double preferredLength);
	PolycrMesh* generateMesh(const CrysesShell& crysesShell, const double preferredLength);
	PolycrMesh* getLastMesh();

	void inputData(string filename);
	void inputData(const CrysesShell& crysesShell);
	void outputData(string filename = "polycr.obj") const;

	Polycrystal3();
	Polycrystal3(string filename);
	Polycrystal3(const CrysesShell& crysesShell);
	~Polycrystal3();
};