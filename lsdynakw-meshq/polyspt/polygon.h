#pragma once
#include "edge.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
struct polygon
{
    static constexpr auto n = 2;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using facet_type = spt::edge<Dim, Real>;

    std::list<facet_type*> edges;

    bool empty() const
    {
        return edges.empty();
    }

    template <typename SubPolytope>
    auto all_of() const
    {
        static_assert(SubPolytope::n == 1 || SubPolytope::n == 0);
        if constexpr (SubPolytope::n == 1)
            return edges;
        else
        {
            std::list<SubPolytope*> res;
            for (const auto& edge : edges)
                for (const auto& vert : edge->vertices)
                    if (std::find(res.begin(), res.end(), vert) == res.end())
                        res.push_back(vert);

            return res;
        }
    }

    bool contains(const facet_type * edge) const
    {
        return edges[0] == edge || edges[1] == edge || edges[2] == edge;
    }

    bool contains(const typename facet_type::facet_type * vert) const
    {
        return edges[0]->contains(vert) || edges[1]->contains(vert) || edges[2]->contains(vert);
    }

    polygon(const polygon& poly)
    {
        edges = poly.edges;
    }
    polygon(const std::list<facet_type*>& edges)
        : edges(edges) {}
    polygon(const facet_type * e0, const facet_type * e1, const facet_type * e2)
        : edges({ e0, e1, e2 }) {}
};

} // namespace spt
