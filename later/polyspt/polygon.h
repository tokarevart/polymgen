#pragma once
#include "edge.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using polygon = polytope<2, Dim, Real>;

} // namespace spt
