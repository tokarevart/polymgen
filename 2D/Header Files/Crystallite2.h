#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Crystallite2
{
public:
	list<ShellEdge*> shellEdges;
	list<Simplex2*> simplexes;
	Polycrystalline* inclInPolycrys;

	// Use when simplexes list already filled.
	bool IsContaining(Simplex2 &simp);
	bool IsContaining(Edge &edge);
	bool IsContaining(Node2 &node);
	// Calculate based on position in space.
	bool Contains(Simplex2 &simp);
	bool Contains(Edge &edge);
	bool Contains(Node2 &node);

	Crystallite2();
	~Crystallite2();
};