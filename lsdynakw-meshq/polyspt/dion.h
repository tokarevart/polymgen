#pragma once
#include <array>
#include "vertex.h"

namespace spt {

template <std::size_t Dim = 3, typename ValueType = typename spt::vec<Dim>::value_type>
using dion = polytope<1, Dim, ValueType>;


template <std::size_t Dim, typename ValueType>
struct polytope<1, Dim, ValueType> {
    static constexpr std::size_t n = 1;
    static constexpr std::size_t dim = Dim;
    using value_type = ValueType;
    using vertex_type = spt::vertex<Dim, ValueType>;
    using facet_type = vertex_type;

    std::array<vertex_type*, 2> vertices = { nullptr, nullptr };

    value_type magnitude() const {
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
        : vertices(vertices) {}
    template <typename... Vertices>
    polytope(Vertices... verts)
        : vertices{ verts... } {}
};

} // namespace spt
