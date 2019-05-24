// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include "data-structs/polyshell.h"


namespace psg {

PolyShell generateCuboids( std::size_t nX, std::size_t nY, std::size_t nZ,
                           real_t dX = 1.0, real_t dY = 1.0, real_t dZ = 1.0 );

} // namespace polygen
