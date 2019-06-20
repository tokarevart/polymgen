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
    using vertex_type = facet_type;

    std::array<vertex_type*, 2> vertices = { nullptr, };

    bool empty() const
    {
        return !vertices[0] && !vertices[1];
    }

    real_type magnitude() const
    {
        return (vertices[1]->pos - vertices[0]->pos).magnitude();
    }

    bool contains(const facet_type* vert) const
    {
        return vertices[0] == vert || vertices[1] == vert;
    }

    polytope(const polytope& poly)
    {
        vertices = poly.vertices;
    }
    polytope(const std::array<facet_type*, 2>& vertices)
        : vertices(vertices) {}
    polytope(const facet_type* v0, const facet_type* v1)
        : vertices({ v0, v1 }) {}
};

} // namespace spt
