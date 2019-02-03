#pragma once
#include <stddef.h>
#include <memory>
#include "../data-structures/PolyStruct.h"


namespace polygen
{
    ::std::unique_ptr<PolyStruct> generateCuboidsPolycrystal(
        size_t nX, size_t nY, size_t nZ,
        double dX = 1.0, double dY = 1.0, double dZ = 1.0);
}
