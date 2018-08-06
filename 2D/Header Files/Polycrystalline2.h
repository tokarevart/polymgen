#pragma once
#include <fstream>
#include <list>
#include <vector>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Polycrystalline2
{
private:
	double minShellNodesCoor[2];
	double maxShellNodesCoor[2];
	double l_min, l_av, l_max;

public:
	list<Crystallite2*> crystallites;
	list<Simplex2*> freeSimplexes;
	vector<Edge2*> freeEdges;
	vector<Node2*> freeNodes;

	void AddCrystallite();                          // May be useless if link from crystallite to polycrystalline will not be used.
	void AddCrystallite(Crystallite2* const& crys); // May be useless if link from crystallite to polycrystalline will not be used.
	size_t GenerateFreeNodesEvenly(double* const& polycrysSizeAxis, size_t* const& minNodesNumAxis);
	void GenerateFreeSimplexesFromFreeNodes();
	void GenerateFreeUniformMesh(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	void FitFreeMeshToShells(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	void DeleteExternalNodes();
	void FillDataWithRemainingDependences(); // 1-st step can be well parallelized by nodes.
	void DistributeNodesEvenly();
	void GenerateFiniteElementMesh(); // Executes other functions.
	// Use when simplexes list already filled.
	const bool IsContaining(const Crystallite2& crys) const;
	const bool IsContaining(const Simplex2& simp) const;
	const bool IsContaining(const Edge2& edge) const;
	const bool IsContaining(const Node2& node) const;
	// Calculate based on position in space.
	const bool Contains(const Simplex2& simp) const;
	const bool Contains(const Edge2& edge) const;
	const bool Contains(const Node2& node) const;
	void InputData(ifstream& shell_nodes, ifstream& crys_edges, ifstream& gener_params);
	void InputData(const vector<double>& shell_nodes, const vector<size_t>& crys_edges, const vector<size_t>& gener_params);
	void OutputData(ofstream& nodesData, ofstream& feNodesData) const;
	void OutputData(vector<double>& nodesData, vector<size_t>& feNodesData) const;

	Polycrystalline2();
	~Polycrystalline2();
};