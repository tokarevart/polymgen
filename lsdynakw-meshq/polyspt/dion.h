#pragma once
#include <array>
#include "vertex.h"

namespace spt {

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::value_type>
using dion = polytope<1, Dim, Real>;

template <typename Real> using dion2 = dion<2, Real>;
using dion2f = dion2<float>;
using dion2d = dion2<double>;
template <typename Real> using dion3 = dion<3, Real>;
using dion3f = dion3<float>;
using dion3d = dion3<double>;

template <std::size_t Dim, typename Real>
struct polytope<1, Dim, Real> {
    static constexpr std::size_t n = 1;
    static constexpr std::size_t dim = Dim;
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, Real>;
    using facet_type = vertex_type;

    std::array<vertex_type*, 2> vertices = { nullptr, nullptr };

    real_type magnitude() const {
        return (vertices[1]->pos - vertices[0]->pos).magnitude();
    }

    bool empty() const {
        return !vertices[0] && !vertices[1];
    }

    bool contains(const vertex_type* vert) const {
        return vertices[0] == vert || vertices[1] == vert;
    }

    polytope(const polytope& poly) {
        vertices = poly.vertices;
    }
    polytope(const std::array<vertex_type*, 2>& vertices)
        : vertices{ vertices } {}
    template <typename... Vertices>
    polytope(Vertices... verts)
        : vertices{ const_cast<vertex_type*>(verts)... } {}
};

} // namespace spt
