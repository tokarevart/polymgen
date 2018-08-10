#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

using std::unique_ptr;
using std::list;


class Crystallite2
{
public:
	list<unique_ptr<ShellEdge2>*> shellEdges;

	// Calculate based on position in space.
	const bool Contains(const Node2& node) const;

	Crystallite2();
	~Crystallite2();
};