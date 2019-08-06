#pragma once
#include "polygon.h"

namespace spt {

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::value_type>
using face = polygon<Dim, Real>;

template <typename Real> using face2 = face<2, Real>;
using face2f = face2<float>;
using face2d = face2<double>;
template <typename Real> using face3 = face<3, Real>;
using face3f = face3<float>;
using face3d = face3<double>;

} // namespace spt
