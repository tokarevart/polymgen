#pragma once
#include "simplex/simplex-2.h"

namespace spt {

template <std::size_t Dim = spt::simplex<2>::dim, typename Real = typename spt::simplex<2, Dim>::real_type>
using triangle = spt::simplex<2, Dim, Real>;

template <std::size_t Dim = spt::simplex<2>::dim, typename Real = typename spt::simplex<2, Dim>::real_type>
using triangle_v = spt::simplex_v<2, Dim, Real>;

} // namespace spt
