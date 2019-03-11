#include "data-structures/polymesh.h"

using pmg::PolyMesh;




PolyMesh::~PolyMesh()
{
    delete[] nodesPositions;
    delete[] tetrs;
    delete[] nCrysesTetrs;
}
