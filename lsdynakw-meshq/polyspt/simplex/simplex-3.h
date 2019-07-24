#pragma once
#include "simplex-2.h"
#include "../algs.h"

namespace spt {

template <std::size_t Dim, typename Real>
struct simplex<3, Dim, Real> {
    static constexpr std::size_t n = 3;
    static constexpr std::size_t dim = Dim;
    using real_type = Real;
    using face_type = spt::simplex<2, Dim, Real>;
    using edge_type = spt::edge<Dim, Real>;
    using vertex_type = spt::vertex<Dim, Real>;
    using facet_type = face_type;

    std::array<face_type*, 4> faces = { nullptr, nullptr, nullptr, nullptr, };

    auto& facets() { return faces; }
    auto& facets() const { return facets(); }

    real_type volume() const {
        auto verts = all_of<vertex_type>();
        return static_cast<real_type>(1.0 / 6.0) *
            std::abs(spt::mixed(
                verts[1]->pos - verts[0]->pos,
                verts[2]->pos - verts[0]->pos,
                verts[3]->pos - verts[0]->pos));
    }

    template <typename SubSimplex>
    auto all_of() const {
        static_assert(
            std::is_same<face_type, SubSimplex>() ||
            std::is_same<edge_type, SubSimplex>() ||
            std::is_same<vertex_type, SubSimplex>());
        if constexpr (std::is_same<face_type, SubSimplex>())
            return faces;
        else if constexpr (std::is_same<edge_type, SubSimplex>()) {
            std::array<edge_type*, 6> edges = { nullptr, };
            std::size_t idx = 0;
            for (const auto& face : faces)
                for (const auto& edge : face->edges)
                    if (std::find(edges.begin(), edges.begin() + idx, edge) == edges.begin() + idx)
                        edges[idx++] = edge;

            return edges;
        } else {
            std::array<vertex_type*, 4> verts = { nullptr, };
            std::size_t idx = 0;
            for (const auto& face : faces)
                for (const auto& edge : face->edges)
                    for (const auto& vert : edge->vertices)
                        if (std::find(verts.begin(), verts.begin() + idx, vert) == verts.begin() + idx)
                            verts[idx++] = vert;

            return verts;
        }
    }

    bool empty() const {
        return !faces[0] && !faces[1] && !faces[2] && !faces[3];
    }

    bool contains(const face_type* face) const {
        return
            faces[0] == face ||
            faces[1] == face ||
            faces[2] == face ||
            faces[3] == face;
    }
    bool contains(const edge_type* edge) const {
        return
            faces[0]->contains(edge) ||
            faces[1]->contains(edge) ||
            faces[2]->contains(edge) ||
            faces[3]->contains(edge);
    }
    bool contains(const vertex_type* vert) const {
        return
            faces[0]->contains(vert) ||
            faces[1]->contains(vert) ||
            faces[2]->contains(vert) ||
            faces[3]->contains(vert);
    }

    simplex(const simplex& simp) {
        faces = simp.faces;
    }
    simplex(const std::array<face_type*, 4>& faces)
        : faces(faces) {}
    template <typename... Faces>
    simplex(Faces... faces)
        : faces{ const_cast<face_type*>(faces)... } {}
};


template <std::size_t Dim, typename Real>
struct simplex_v<3, Dim, Real> {
    static constexpr std::size_t n = 3;
    static constexpr std::size_t dim = Dim;
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, Real>;

    std::array<vertex_type*, 4> vertices = { nullptr, nullptr, nullptr, nullptr };

    real_type volume() const {
        return static_cast<real_type>(1.0 / 6.0) *
            std::abs(spt::mixed(
                vertices[1]->pos - vertices[0]->pos,
                vertices[2]->pos - vertices[0]->pos,
                vertices[3]->pos - vertices[0]->pos));
    }

    bool empty() const {
        return !vertices[0] && !vertices[1] && !vertices[2] && !vertices[3];
    }

    bool contains(const vertex_type* vert) const {
        return std::find(vertices.begin(), vertices.end(), vert) != vertices.end();
    }

    simplex_v(const simplex_v& other) {
        vertices = other.vertices;
    }
    simplex_v(const std::array<vertex_type*, 4> & verts)
        : vertices(vertices) {}
    template <typename... Vertices>
    simplex_v(Vertices... verts)
        : vertices{ const_cast<vertex_type*>(verts)... } {}
};

} // namespace spt
