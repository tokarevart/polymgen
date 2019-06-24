#pragma once
#include "simplex/simplex-2.h"

namespace spt {

template <std::size_t Dim, typename Real>
using triangle = spt::simplex<2, Dim, Real>;

template <std::size_t Dim, typename Real>
using triangle_v = spt::simplex_v<2, Dim, Real>;

} // namespace spt
