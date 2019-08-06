#pragma once
#include "simplex/simplex-3.h"

namespace spt {

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::value_type>
using tetrahedron = spt::simplex<3, Dim, Real>;

template <typename Real> using tetrahedron2 = tetrahedron<2, Real>;
using tetrahedron2f = tetrahedron2<float>;
using tetrahedron2d = tetrahedron2<double>;
template <typename Real> using tetrahedron3 = tetrahedron<3, Real>;
using tetrahedron3f = tetrahedron3<float>;
using tetrahedron3d = tetrahedron3<double>;

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::value_type>
using tetrahedron_v = spt::simplex_v<3, Dim, Real>;

template <typename Real> using tetrahedron2_v = tetrahedron_v<2, Real>;
using tetrahedron2f_v = tetrahedron2_v<float>;
using tetrahedron2d_v = tetrahedron2_v<double>;
template <typename Real> using tetrahedron3_v = tetrahedron_v<3, Real>;
using tetrahedron3f_v = tetrahedron3_v<float>;
using tetrahedron3d_v = tetrahedron3_v<double>;

} // namespace spt
