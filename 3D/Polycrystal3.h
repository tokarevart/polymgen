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

	ShellEdge3* FindShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const;
	unique_ptr<Edge3>* FindStartFrontEdge(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1) const;

	template <class T>
	void ErasePtrsToNullptr(vector<unique_ptr<T>*>& vec);

	void SetLinksWithShell();
	void ErasePtrsToNullptrFromVectors();
	void TriangulateShell();
	void StartFrontDelaunayPostprocessing();

	void InputData();

public:
	void TriangulatePolycrystalNoStruct(const double preferredLength);
	void TriangulatePolycrystalNoStruct(string filename, const double preferredLength);
	void TriangulatePolycrystalNoStruct(const CrysesShell& crysesShell, const double preferredLength);
	PolycrMesh* StructurizeTriangulation();
	PolycrMesh* TriangulatePolycrystal(const double preferredLength);
	PolycrMesh* TriangulatePolycrystal(string filename, const double preferredLength);
	PolycrMesh* TriangulatePolycrystal(const CrysesShell& crysesShell, const double preferredLength);
	PolycrMesh* GetLastTriangulation();

	void InputData(string filename);
	void InputData(const CrysesShell& crysesShell);
	void OutputData(string filename = "polycr.obj") const;

	Polycrystal3();
	Polycrystal3(string filename);
	Polycrystal3(const CrysesShell& crysesShell);
	~Polycrystal3();
};