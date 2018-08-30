#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

using std::unique_ptr;
using std::list;

class ShellEdge3
{
public:
	ShellVertex3* vertexes[2];
	vector<ShellFacet3*> inclInFacets;

	const double Magnitude() const;
	const double SqrMagnitude() const;

	const bool IsContaining(const ShellVertex3& vertex) const;

	ShellEdge3();
	ShellEdge3(ShellVertex3& vertex0, ShellVertex3& vertex1);
	~ShellEdge3();
};