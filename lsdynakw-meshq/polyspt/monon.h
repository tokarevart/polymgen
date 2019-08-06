#pragma once
#include "polytope-base.h"

namespace spt {

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::value_type>
using monon = polytope<0, Dim, Real>;

template <typename Real> using monon2 = monon<2, Real>;
using monon2f = monon2<float>;
using monon2d = monon2<double>;
template <typename Real> using monon3 = monon<3, Real>;
using monon3f = monon3<float>;
using monon3d = monon3<double>;

template <std::size_t Dim, typename Real>
struct polytope<0, Dim, Real> {
    static constexpr std::size_t n = 0;
    static constexpr std::size_t dim = Dim;
    using value_type = Real;

    spt::vec<Dim, Real> pos;

    polytope(const polytope& other)
        : pos{ other.pos } {}
    polytope(const spt::vec<Dim, Real>& pos)
        : pos{ pos } {}
};

} // namespace spt
