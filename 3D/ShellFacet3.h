#pragma once
#include <vector>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

class ShellFacet3
{
public:
	ShellEdge3* edges[3];

	const bool IsContaining(const ShellEdge3& edge) const;
	const bool IsContaining(const ShellVertex3& vertex) const;

	ShellFacet3();
	ShellFacet3(
		ShellEdge3& edge0, 
		ShellEdge3& edge1, 
		ShellEdge3& edge2);
	~ShellFacet3();
};