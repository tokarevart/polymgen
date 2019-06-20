#pragma once
#include "vertex.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using dion = polytope<1, Dim, Real>;


template <std::size_t Dim, typename Real>
struct polytope<1, Dim, Real>
{
    static constexpr auto n = 1;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using facet_type = spt::vertex<Dim, Real>;

    std::array<facet_type*, 2> facets = { nullptr, };

    bool empty() const
    {
        return !facets[0] && !facets[1];
    }

    template <typename SubPolytope>
    std::array<SubPolytope*, 2> all_of() const
    {
        static_assert(SubPolytope::n == 0);
        return facets;
    }

    bool contains(const facet_type* subpt) const
    {
        return facets[0] == subpt || facets[1] == subpt;
    }

    polytope(const polytope& poly)
    {
        facets = poly.facets;
    }
    polytope(const std::array<facet_type*, 2>& facets)
        : facets(facets) {}
    polytope(const facet_type* f0, const facet_type* f1)
        : facets({ f0, f1 }) {}
};

} // namespace spt
