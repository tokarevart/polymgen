#pragma once
#include "monon.h"

namespace spt {

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::value_type>
using vertex = monon<Dim, Real>;

template <typename Real> using vertex2 = vertex<2, Real>;
using vertex2f = vertex2<float>;
using vertex2d = vertex2<double>;
template <typename Real> using vertex3 = vertex<3, Real>;
using vertex3f = vertex3<float>;
using vertex3d = vertex3<double>;

} // namespace spt
