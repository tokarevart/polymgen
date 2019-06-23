#pragma once
#include "polyhedron.h"

namespace spt {

template <std::size_t Dim = 3, typename ValueType = typename spt::vec<Dim>::value_type>
using cell = polyhedron<Dim, ValueType>;

} // namespace spt
