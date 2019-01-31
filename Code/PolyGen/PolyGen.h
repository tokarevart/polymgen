#pragma once
#include "../Data Structures/PolyStruct.h"


namespace polygen
{
    PolyStruct* generateCuboidsPolycrystal(
        size_t nX, size_t nY, size_t nZ, 
        double dX = 1.0, double dY = 1.0, double dZ = 1.0);
}