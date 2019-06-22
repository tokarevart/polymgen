#pragma once
#include "edge.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using polygon = polytope<2, Dim, Real>;

template <std::size_t Dim, typename Real>
struct polytope<2, Dim, Real>
{
    static constexpr auto n = 2;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, Real>;
    using edge_type = spt::edge<Dim, Real>;
    using facet_type = edge_type;

    std::list<edge_type*> edges;

    bool empty() const
    {
        return edges.empty();
    }

    template <typename SubPolytope>
    auto all_of() const
    {
        static_assert(
            std::is_same<edge_type, SubPolytope>() ||
            std::is_same<vertex_type, SubPolytope>());
        if constexpr (std::is_same<edge_type, SubPolytope>())
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

    bool contains(const edge_type* edge) const
    {
        return edges[0] == edge || edges[1] == edge || edges[2] == edge;
    }

    bool contains(const vertex_type* vert) const
    {
        return edges[0]->contains(vert) || edges[1]->contains(vert) || edges[2]->contains(vert);
    }

    polytope(const polytope& poly)
    {
        edges = poly.edges;
    }
    polytope(const std::list<edge_type*>& edges)
        : edges(edges) {}
    template <typename... Edges>
    polytope(Edges... edges)
        : edges{ edges... } {}
};

} // namespace spt
