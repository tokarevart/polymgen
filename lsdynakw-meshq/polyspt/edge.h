#pragma once
#include "dion.h"

namespace spt {

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::value_type>
using edge = dion<Dim, Real>;

template <typename Real> using edge2 = edge<2, Real>;
using edge2f = edge2<float>;
using edge2d = edge2<double>;
template <typename Real> using edge3 = edge<3, Real>;
using edge3f = edge3<float>;
using edge3d = edge3<double>;

} // namespace spt
