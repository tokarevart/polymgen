#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

using std::vector;


class Crystallite3
{
public:
	vector<ShellFacet3*> shellFacets;
	vector<ShellEdge3*> shellEdges;

	// Calculate based on position in space.
	const bool Contains(const Vertex3& vertex) const;

	Crystallite3();
	~Crystallite3();
};