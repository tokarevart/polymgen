// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include "../data-structs/polyshell.h"


namespace psg {

PolyShell cuboids(std::size_t nX, std::size_t nY, std::size_t nZ,
                  real_t dX = static_cast<real_t>(1), 
                  real_t dY = static_cast<real_t>(1), 
                  real_t dZ = static_cast<real_t>(1));

} // namespace polygen
