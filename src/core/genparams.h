#pragma once
#include <cstddef>
#include "../real-type.h"

namespace pmg {
namespace genparams {

struct Surface {
    real_t min_dis = static_cast<real_t>(0.2);
    real_t max_def = static_cast<real_t>(0.3);
    real_t def = static_cast<real_t>(0.5);

    std::size_t n_smooth_iters = 20;
    std::size_t n_delaunay_smooth_iters = 3;
};

using Shell = Surface;

struct Volume {
    real_t min_dis = static_cast<real_t>(0.2);
    real_t max_def = static_cast<real_t>(0.3);
    real_t def = static_cast<real_t>(0.4);

    std::size_t n_smooth_iters = 20;
};

struct Polyhedron {
    Shell  shell = Shell();
    Volume volume = Volume();
};

} // namespace genparams
} // namespace pmg
