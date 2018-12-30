#pragma once
#include <list>
#include <memory>
#include "Definitions.h"
#include "Inclusions.h"

using std::unique_ptr;
using std::list;

class ShellEdge3
{
public:
	ShellVertex3* vertexes[2];

	const double magnitude() const;
	const double sqrMagnitude() const;

	const bool contains(const ShellVertex3& vertex) const;

	ShellEdge3();
	ShellEdge3(ShellVertex3& vertex0, ShellVertex3& vertex1);
	~ShellEdge3();
};