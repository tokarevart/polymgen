#pragma once
#include "vertex.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using edge = polytope<1, Dim, Real>;

template <std::size_t Dim, typename Real>
class polytope<1, Dim, Real>
{
public:
    static constexpr auto n = 1;
    static constexpr auto dim = Dim;
    using real_type = Real;

    const std::array<vertex<Dim, Real>*, 2> facets;

    template <typename SubPolytope>
    std::array<SubPolytope*, 2> all_of() const
    {
        if constexpr (SubPolytope::n >= 1)
            throw std::logic_error("SubPolytope::n must be less than this class n");
        else
            return facets;
    }
    
    bool contains(const vertex<Dim, Real>* subpt) const
    {
        return facets[0] == subpt || facets[1] == subpt;
    }

    polytope() = delete;
    polytope(const std::array<vertex<Dim, Real>*, 2>& facets)
        : facets(facets) {}
    polytope(const vertex<Dim, Real>* v0, const vertex<Dim, Real>* v1, const vertex<Dim, Real>* v2)
        : facets({ v0, v1, v2 }) {}
};

} // namespace spt
