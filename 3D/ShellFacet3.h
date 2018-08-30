#pragma once
#include <vector>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

class ShellFacet3
{
public:
	ShellEdge3* edges[3];

	vector<unique_ptr<Edge3>*> innerEdges;

	ShellFacet3();
	ShellFacet3(ShellEdge3& edge0, ShellEdge3& edge1, ShellEdge3& edge2);
	~ShellFacet3();
};