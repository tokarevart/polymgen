#pragma once
#include <fstream>
#include <list>
#include <vector>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Polycrystalline
{
public:
	list<Crystallite2*> crystallites;

	void AddCrystallite();                    // May be useless.
	void AddCrystallite(Crystallite2* &crys); // May be useless.
	// Use when simplexes list already filled.
	bool IsContaining(Simplex2 &simp);
	bool IsContaining(Edge &edge);
	bool IsContaining(Node2 &node);
	// Calculate based on position in space.
	bool Contains(Crystallite2 &crys);
	bool Contains(Simplex2 &simp);
	bool Contains(Edge &edge);
	bool Contains(Node2 &node);
	void InputData(ifstream &shell_nodes, ifstream &crys_planes);
	void InputData(vector<double> &shell_nodes, vector<size_t> &crys_planes);
	void OutputData(ofstream &nodesData, ofstream &feNodesData);
	void OutputData(vector<double> &nodesData, vector<size_t> &feNodesData);

	Polycrystalline();
	~Polycrystalline();
};