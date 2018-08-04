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
	double l_min, l_av, l_max;  // l_av must be less than the least shell edge.

public:
	list<Crystallite2*> crystallites;
	list<Simplex2*> freeSimplexes;
	list<Node2*> freeNodes;

	void AddCrystallite();                          // May be useless if link from crystallite to polycrystalline will not be used.
	void AddCrystallite(Crystallite2* const& crys); // May be useless if link from crystallite to polycrystalline will not be used.
	void GenerateFreeUniformMesh(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	void FitFreeMeshToShells(); // There may will not be significant performance improvement. But if you want you can parallize 1-st step.
	void DeleteExternalNodes();
	void FillDataWithRemainingDependences(); // 1-st step can be well parallelized by nodes.
	void DistributeNodesEvenly();
	void GenerateFiniteElementMesh(); // Executes other functions.
	// Use when simplexes list already filled.
	const bool IsContaining(const Crystallite2& crys);
	const bool IsContaining(const Simplex2& simp);
	const bool IsContaining(const Edge2& edge);
	const bool IsContaining(const Node2& node);
	// Calculate based on position in space.
	const bool Contains(const Simplex2& simp);
	const bool Contains(const Edge2& edge);
	const bool Contains(const Node2& node);
	void InputData(ifstream& shell_nodes, ifstream& crys_edges, ifstream& gener_params);
	void InputData(const vector<double>& shell_nodes, const vector<size_t>& crys_edges, const vector<size_t>& gener_params);
	void OutputData(ofstream& nodesData, ofstream& feNodesData);
	void OutputData(vector<double>& nodesData, vector<size_t>& feNodesData);

	Polycrystalline2();
	~Polycrystalline2();
};