#pragma once
#include "polytope-base.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using monon = polytope<0, Dim, Real>;


template <std::size_t Dim, typename Real>
struct polytope<0, Dim, Real> {
    static constexpr std::size_t n = 0;
    static constexpr std::size_t dim = Dim;
    using real_type = Real;

    spt::vec<Dim, Real> pos;

    polytope(const polytope& other)
        : pos(other.pos) {}
    polytope(const spt::vec<Dim, Real>& pos)
        : pos(pos) {}
};

} // namespace spt
