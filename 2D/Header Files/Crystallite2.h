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
	const bool IsContaining(const Simplex2& simp) const;
	const bool IsContaining(const Edge2& edge) const;
	const bool IsContaining(const Node2& node) const;
	// Calculate based on position in space.
	const bool Contains(const Simplex2& simp) const;
	const bool Contains(const Edge2& edge) const;
	const bool Contains(const Node2& node) const;

	Crystallite2();
	~Crystallite2();
};