#include "CrysesShell.h"

CrysesShell::CrysesShell() {}

CrysesShell::~CrysesShell()
{
    delete[] nodesPositions;
    delete[] facets;
    delete[] nCrysesFacets;
    delete[] cryses;
}