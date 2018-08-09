#pragma once
#include <list>
#include <memory>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Crystallite2
{
public:
	std::list<std::unique_ptr<ShellEdge2>*> shellEdges;
	std::list<std::unique_ptr<Simplex2>*> simplexes;

	// Calculate based on position in space.
	const bool Contains(const Node2& node) const;

	Crystallite2();
	~Crystallite2();
};