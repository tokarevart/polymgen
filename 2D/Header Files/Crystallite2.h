#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Crystallite2
{
public:
	list<ShellEdge2*> shellEdges;
	list<Simplex2*> simplexes;
	Polycrystalline2* inclInPolycrys;

	// Use when simplexes list already filled.
	bool IsContaining(Simplex2 &simp);
	bool IsContaining(Edge2 &edge);
	bool IsContaining(Node2 &node);
	// Calculate based on position in space.
	bool Contains(Simplex2 &simp);
	bool Contains(Edge2 &edge);
	bool Contains(Node2 &node);

	Crystallite2();
	~Crystallite2();
};