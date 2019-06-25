#pragma once
#include "polygon.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::value_type>
using face = polygon<Dim, Real>;

} // namespace spt
