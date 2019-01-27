#include "PolyMesh.h"

PolyMesh::PolyMesh() {}

PolyMesh::~PolyMesh()
{
    delete[] nodesPositions;
    delete[] tetrs;
    delete[] nCrysesTetrs;
}