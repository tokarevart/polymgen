#pragma once
#include "simplex/simplex-2.h"

namespace spt {

template <std::size_t Dim = spt::simplex<2>::dim, typename Real = typename spt::simplex<2, Dim>::real_type>
using triangle = spt::simplex<2, Dim, Real>;

template <typename Real> using triangle2 = triangle<2, Real>;
using triangle2f = triangle2<float>;
using triangle2d = triangle2<double>;
template <typename Real> using triangle3 = triangle<3, Real>;
using triangle3f = triangle3<float>;
using triangle3d = triangle3<double>;

template <std::size_t Dim = spt::simplex<2>::dim, typename Real = typename spt::simplex<2, Dim>::real_type>
using triangle_v = spt::simplex_v<2, Dim, Real>;

template <typename Real> using triangle2_v = triangle_v<2, Real>;
using triangle2f_v = triangle2_v<float>;
using triangle2d_v = triangle2_v<double>;
template <typename Real> using triangle3_v = triangle_v<3, Real>;
using triangle3f_v = triangle3_v<float>;
using triangle3d_v = triangle3_v<double>;

} // namespace spt
