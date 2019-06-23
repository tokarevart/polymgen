#pragma once
#include <cstddef>
#include "../real-type.h"

namespace pmg {
namespace genparams {

struct Surface {
    real_t cMinDis = static_cast<real_t>(0.2);
    real_t cMaxDef = static_cast<real_t>(0.3);
    real_t cDef = static_cast<real_t>(0.5);

    std::size_t nSmoothIters = 20;
    std::size_t nDelaunaySmoothIters = 3;
};

using Shell = Surface;

struct Volume {
    real_t cMinDis = static_cast<real_t>(0.2);
    real_t cMaxDef = static_cast<real_t>(0.3);
    real_t cDef = static_cast<real_t>(0.4);

    std::size_t nSmoothIters = 20;
};

struct Polyhedron {
    Shell  shell = Shell();
    Volume volume = Volume();
};

} // namespace genparams
} // namespace pmg
