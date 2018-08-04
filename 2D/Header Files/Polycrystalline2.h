#pragma once
#include <fstream>
#include <list>
#include <vector>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Polycrystalline2
{
public:
	list<Crystallite2*> crystallites;
	list<Simplex2*> freeSimplexes; // In GenerateUniformGrid() also use local list<Node2*> or vector<Node2*>

	void AddCrystallite();                    // May be useless.
	void AddCrystallite(Crystallite2* &crys); // May be useless.
	// Use when simplexes list already filled.
	bool IsContaining(Crystallite2 &crys);
	bool IsContaining(Simplex2 &simp);
	bool IsContaining(Edge2 &edge);
	bool IsContaining(Node2 &node);
	// Calculate based on position in space.
	bool Contains(Simplex2 &simp);
	bool Contains(Edge2 &edge);
	bool Contains(Node2 &node);
	void InputData(ifstream &shell_nodes, ifstream &crys_edges);
	void InputData(vector<double> &shell_nodes, vector<size_t> &crys_edges);
	void OutputData(ofstream &nodesData, ofstream &feNodesData);
	void OutputData(vector<double> &nodesData, vector<size_t> &feNodesData);

	Polycrystalline2();
	~Polycrystalline2();
};