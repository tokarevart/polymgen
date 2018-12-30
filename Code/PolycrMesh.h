#pragma once

struct PolycrMesh
{
	size_t  nodesNum;

	// { x0, y0, z0, x1, y1, z1 ... }
	double* nodesPositions;
	size_t  tetrsNum;
	size_t  crysesNum;

	// Number of tetrahedrons in each crystalline.
	size_t* crysesTetrsNum;

	// c - crystallite,
	// t - tetrahedron,
	// n - node.
	// { ... c[i]_t[0]_n[0..3] ... c[i]_t[crysesTetrsNum[i]]_n[0..3] }
	size_t* tetrs;

	PolycrMesh();
	~PolycrMesh();
};