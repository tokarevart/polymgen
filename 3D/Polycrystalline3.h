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
	double _preferredLength;

	vector<ShellFacet3*> _shellFacets;
	vector<ShellEdge3*> _shellEdges;
	vector<ShellVertex3*> _shellVertexes;

	vector<unique_ptr<Facet3>*> _startFrontFacets;
	vector<unique_ptr<Edge3>*> _startFrontEdges;
	vector<unique_ptr<Vertex3>*> _startFrontVertexes;

	template <class T>
	void ErasePtrsToNullptr(vector<unique_ptr<T>*> &vec);

public:
	vector<Crystallite3*> crystallites;

	//void SetShellFacetsInnerEdges();
	void SetLinksWithShell();
	//unique_ptr<Edge3>* FindEdge(Vertex3& vertex0, Vertex3& vertex1);
	void ErasePtrsToNullptrFromVectors();
	void TriangulateShell();
	//Vector3 ShiftToLavEdges(const Vertex3& vertex);
	//void DistributeVertexesOnShellEvenly(size_t iterations_num = 0);
	void DelaunayPostprocessing();
	void TriangulateCrystallitesVolumes();
	//void GenerateFiniteElementMesh(); // Executes other functions.

	void InputData();
	
	void OutputData(string filename = "polycr.obj") const;

	Polycrystalline3();
	~Polycrystalline3();
};