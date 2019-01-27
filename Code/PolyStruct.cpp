#include "PolyStruct.h"

PolyStruct::PolyStruct() {}

PolyStruct::~PolyStruct()
{
    delete[] nodesPositions;
    delete[] facets;
    delete[] nCrysesFacets;
    delete[] cryses;
}