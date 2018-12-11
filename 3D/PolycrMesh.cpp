#include "PolycrMesh.h"

PolycrMesh::PolycrMesh() {}

PolycrMesh::~PolycrMesh()
{
	delete[] nodesPositions;
	delete[] tetrs;
	delete[] crysesTetrsNum;
}