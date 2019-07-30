#pragma once
#include "simplex-base.h"

namespace spt {

template <std::size_t Dim, typename Real>
struct simplex<2, Dim, Real> {
    static constexpr std::size_t n = 2;
    static constexpr std::size_t dim = Dim;
    using real_type = Real;
    using edge_type = spt::edge<Dim, Real>;
    using vertex_type = spt::vertex<Dim, Real>;
    using facet_type = edge_type;

    std::array<edge_type*, 3> edges{ nullptr, nullptr, nullptr };

    auto& facets() { return edges; }
    auto& facets() const { return facets(); }

    template <typename SubSimplex>
    auto all_of() const {
        static_assert(std::is_same<edge_type, SubSimplex>() || std::is_same<vertex_type, SubSimplex>());
        if constexpr (std::is_same<edge_type, SubSimplex>())
            return edges;
        else {
            std::array<SubSimplex*, 4> verts = { nullptr, };
            std::size_t idx = 0;
            for (const auto& edge : edges)
                for (const auto& vert : edge->vertices)
                    if (std::find(verts.begin(), verts.begin() + idx, vert) == verts.begin() + idx)
                        verts[idx++] = vert;

            return verts;
        }
    }

    real_type area() const {
        vec<3> vec0 = edges[0]->vertices[1]->pos - edges[0]->vertices[0]->pos;
        vec<3> vec1 = edges[1]->vertices[1]->pos - edges[1]->vertices[0]->pos;
        return static_cast<real_type>(0.5) * spt::cross(vec0, vec1).magnitude();
    }

    bool empty() const {
        return !edges[0] && !edges[1] && !edges[2];
    }

    bool contains(const edge_type* edge) const {
        return edges[0] == edge || edges[1] == edge || edges[2] == edge;
    }
    bool contains(const vertex_type* vert) const {
        return edges[0]->contains(vert) || edges[1]->contains(vert) || edges[2]->contains(vert);
    }

    simplex(const simplex& simp) {
        edges = simp.edges;
    }
    simplex(const std::array<edge_type*, 3>& edges)
        : edges{ edges } {}
    template <typename... Edges>
    simplex(Edges... edges)
        : edges{ const_cast<edge_type*>(edges)... } {}
};


template <std::size_t Dim, typename Real>
struct simplex_v<2, Dim, Real> {
    static constexpr std::size_t n = 2;
    static constexpr std::size_t dim = Dim;
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, Real>;

    std::array<vertex_type*, 3> vertices{ nullptr, nullptr, nullptr };

    real_type area() const {
        vec<dim> vec0 = vertices[1]->pos - vertices[0]->pos;
        vec<dim> vec1 = vertices[2]->pos - vertices[0]->pos;

        auto cross = vec<dim>::cross(vec0, vec1);
        real_type doubled_area;
        if constexpr (dim == 2)
            doubled_area = cross;
        else
            doubled_area = cross.magnitude();

        return doubled_area * static_cast<real_type>(0.5);
    }

    bool empty() const {
        return !vertices[0] && !vertices[1] && !vertices[2];
    }

    bool contains(const vertex_type* vert) const {
        return std::find(vertices.begin(), vertices.end(), vert) != vertices.end();
    }

    simplex_v(const simplex_v& other) {
        vertices = other.vertices;
    }
    simplex_v(const std::array<vertex_type*, 3>& verts)
        : vertices{ vertices } {}
    template <typename... Vertices>
    simplex_v(Vertices... verts)
        : vertices{ const_cast<vertex_type*>(verts)... } {}
};

} // namespace spt
