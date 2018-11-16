#pragma once

struct CrystallitesShell
{
	size_t  nodesNum;
	double* nodesPositions;   // { x0, y0, z0, x1, y1, z1 ... }
	size_t  facetsNum;
	size_t* facets;           // { facet0_node0, facet0_node1, facet0_node2, facet1_node0, facet1_node1, facet1_node2 ... }
	size_t  crysesNum;
	size_t* crysesFacetsNums; // Number of facets in each crystallite shell.
	size_t* cryses;           // { crys0_facet0, crys0_facet1, ... crys1_facet0, crys1_facet1 ... }

	CrystallitesShell();
	~CrystallitesShell();
};