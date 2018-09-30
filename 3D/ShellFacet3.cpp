#include "ShellFacet3.h"


const bool ShellFacet3::IsContaining(const ShellEdge3& edge) const
{
	for (auto &edge_ : edges)
	{
		if (edge_ == &edge)
		{
			return true;
		}
	}

	return false;
}

const bool ShellFacet3::IsContaining(const ShellVertex3& vertex) const
{
	for (auto &edge : edges)
	{
		if (edge->IsContaining(vertex))
		{
			return true;
		}
	}

	return false;
}

ShellFacet3::ShellFacet3() {}

ShellFacet3::ShellFacet3(ShellEdge3& edge0, ShellEdge3& edge1, ShellEdge3& edge2)
{
	edges[0] = &edge0;
	edges[1] = &edge1;
	edges[2] = &edge2;
}

ShellFacet3::~ShellFacet3() {}