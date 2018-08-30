#include "ShellFacet3.h"


ShellFacet3::ShellFacet3()
{
}

ShellFacet3::ShellFacet3(ShellEdge3& edge0, ShellEdge3& edge1, ShellEdge3& edge2)
{
	edges[0] = &edge0;
	edges[1] = &edge1;
	edges[2] = &edge2;
}

ShellFacet3::~ShellFacet3()
{
}