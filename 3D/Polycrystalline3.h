#pragma once
#include <string>
#include <fstream>
#include <vector>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

using std::ifstream;
using std::ofstream;
using std::unique_ptr;
using std::vector;
using std::string;

class Polycrystalline3
{
private:
	double preferredLength;

	// Add external polycrystalline shell and some method for its parsing from inclInCrysesNum
	vector<ShellFacet3*> _shellFacets;
	vector<ShellEdge3*> _shellEdges;
	vector<ShellVertex3*> _shellVertexes;

	vector<unique_ptr<Facet3>*> _freeFacets;
	vector<unique_ptr<Edge3>*> _freeEdges;
	vector<unique_ptr<Vertex3>*> _freeVertexes;

public:
	vector<Crystallite3*> crystallites;

	//void SetShellFacetsInnerEdges();
	void SetLinksWithShell();
	unique_ptr<Edge3>* FindEdge(Vertex3& vertex0, Vertex3& vertex1);
	void ErasePtrsToNullptrFromVectors();
	void TriangulateShell();
	Vector3 ShiftToLavEdges(const Vertex3& vertex);
	void DistributeVertexesOnShellEvenly(size_t iterations_num = 0);
	void DelaunayPostprocessing();
	//void GenerateFiniteElementMesh(); // Executes other functions.

	void InputData();
	
	void OutputData(string filename = "polycr.obj") const;

	Polycrystalline3();
	~Polycrystalline3();
};