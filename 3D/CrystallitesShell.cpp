#include "CrystallitesShell.h"

CrystallitesShell::CrystallitesShell() {}

CrystallitesShell::~CrystallitesShell()
{
	delete[] nodesPositions;
	delete[] facets;
	delete[] crysesFacetsNums;
	delete[] cryses;
}