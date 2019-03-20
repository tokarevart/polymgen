// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <stddef.h>
#include <memory>
#include "data-structures/polystruct.h"


namespace polygen {

PolyStruct generateCuboidsPolyhedralSet(
    size_t nX, size_t nY, size_t nZ,
    real_t dX = 1.0, real_t dY = 1.0, real_t dZ = 1.0) noexcept;

} // namespace polygen
