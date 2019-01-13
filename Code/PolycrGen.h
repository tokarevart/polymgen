#pragma once
#include "CrysesShell.h"


namespace polycrgen
{
    CrysesShell* generateCuboidsPolycrystal(
        size_t nX, size_t nY, size_t nZ, 
        double dX = 1.0, double dY = 1.0, double dZ = 1.0);
}
