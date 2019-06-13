#pragma once
#include "vertex.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using edge = polytope<1, Dim, Real>;


template <std::size_t Dim, typename Real>
struct polytope<1, Dim, Real>
{
    static constexpr auto n = 1;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using facet_type = spt::vertex<Dim, Real>;

    const std::array<facet_type*, 2> facets;

    template <typename SubPolytope>
    std::array<SubPolytope*, 2> all_of() const
    {
        if constexpr (SubPolytope::n >= 1)
            throw std::logic_error("SubPolytope::n must be less than this class n");
        else
            return facets;
    }
    
    bool contains(const facet_type* subpt) const
    {
        return facets[0] == subpt || facets[1] == subpt;
    }

    polytope() = delete;
    polytope(const std::array<facet_type*, 2>& facets)
        : facets(facets) {}
    polytope(const facet_type* v0, const facet_type* v1, const facet_type* v2)
        : facets({ v0, v1, v2 }) {}
};

} // namespace spt
