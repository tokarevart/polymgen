#include "PolycrystalTriangulation.h"

PolycrystalTriangulation::PolycrystalTriangulation() {}

PolycrystalTriangulation::~PolycrystalTriangulation()
{
	delete[] nodesPositions;
	delete[] tetrs;
	delete[] crysesTetrsNum;
}