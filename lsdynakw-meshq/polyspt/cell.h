#pragma once
#include "polyhedron.h"

namespace spt {

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::value_type>
using cell = polyhedron<Dim, Real>;

template <typename Real> using cell2 = cell<2, Real>;
using cell2f = cell2<float>;
using cell2d = cell2<double>;
template <typename Real> using cell3 = cell<3, Real>;
using cell3f = cell3<float>;
using cell3d = cell3<double>;

} // namespace spt
