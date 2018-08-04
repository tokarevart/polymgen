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
	double minShellNodesCoor[2]; // Add its calculation
	double maxShellNodesCoor[2]; // Add its calculation
	double l_min, l_av, l_max;   // Add its calculation

public:
	list<Crystallite2*> crystallites;
	list<Simplex2*> freeSimplexes; // In GenerateUniformGrid() also use local list<Node2*> or vector<Node2*>

	void AddCrystallite();                          // May be useless if link from crystallite to polycrystalline will not be used.
	void AddCrystallite(Crystallite2* const& crys); // May be useless if link from crystallite to polycrystalline will not be used.
	void GenerateUniformGrid();
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