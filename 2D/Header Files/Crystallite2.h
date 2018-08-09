#pragma once
#include <list>
#include <memory>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Crystallite2
{
public:
	std::list<std::unique_ptr<ShellEdge2>*> shellEdges;
	std::list<Simplex2*> simplexes;

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