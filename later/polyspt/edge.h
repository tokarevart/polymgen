#pragma once
#include "dion.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using edge = dion<Dim, Real>;

} // namespace spt
