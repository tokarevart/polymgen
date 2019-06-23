#pragma once
#include "polytope-base.h"

namespace spt {

template <std::size_t Dim = 3, typename ValueType = typename spt::vec<Dim>::value_type>
using monon = polytope<0, Dim, ValueType>;


template <std::size_t Dim, typename ValueType>
struct polytope<0, Dim, ValueType> {
    static constexpr std::size_t n = 0;
    static constexpr std::size_t dim = Dim;
    using value_type = ValueType;

    spt::vec<Dim, ValueType> pos;

    polytope(const polytope& other)
        : pos(other.pos) {}
    polytope(const spt::vec<Dim, ValueType>& pos)
        : pos(pos) {}
};

} // namespace spt
