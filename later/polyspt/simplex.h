#pragma once
#include <cstddef>
#include <array>
#include <algorithm>
#include "face.h"

namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
struct simplex;

template <std::size_t Dim, typename Real>
struct simplex<3, Dim, Real>
{
    static constexpr auto n = 3;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using facet_type = spt::face<Dim, Real>;
    using face_type = facet_type;
    using edge_type = typename face_type::edge_type;
    using vertex_type = typename edge_type::vertex_type;

    std::array<face_type*, 4> faces = { nullptr, };

    bool empty() const
    {
        return !faces[0] && !faces[1] && !faces[2] && !faces[3];
    }

    template <typename SubPolytope>
    auto all_of() const
    {
        static_assert(
            std::is_same<face_type, SubPolytope>() ||
            std::is_same<edge_type, SubPolytope>() ||
            std::is_same<vertex_type, SubPolytope>());
        if constexpr (std::is_same<face_type, SubPolytope>())
            return faces;
        else if constexpr (std::is_same<edge_type, SubPolytope>())
        {
            std::array<edge_type*, 6> edges = { nullptr, };
            std::size_t idx = 0;
            for (const auto& face : faces)
                for (const auto& edge : face->edges)
                    if (std::find(edges.begin(), edges.begin() + idx, edge) == edges.begin() + idx)
                        edges[idx++] = edge;

            return edges;
        }
        else
        {
            std::array<vertex_type*, 4> verts = { nullptr, };
            std::size_t idx = 0;
            for (const auto& face : faces)
                for (const auto& edge : face->edges)
                {
                    if (idx == 4) 
                        return verts;
                    else
                    {
                        for (const auto& vert : edge->vertices)
                            if (std::find(verts.begin(), verts.begin() + idx, vert) == verts.begin() + idx)
                                verts[idx++] = vert;
                    }
                }
                    
            return verts;
        }
    }

    real_type volume() const
    {
        auto verts = all_of<vertex_type>();
        return static_cast<real_type>(1.0 / 6.0) * 
            std::abs(vec<3>::mixed(
                verts[1]->pos - verts[0]->pos, 
                verts[2]->pos - verts[0]->pos, 
                verts[3]->pos - verts[0]->pos));
    }

    bool contains(const face_type* face) const
    {
        return 
            faces[0] == face || 
            faces[1] == face || 
            faces[2] == face || 
            faces[3] == face;
    }

    bool contains(const edge_type* edge) const
    {
        return 
            faces[0]->contains(edge) || 
            faces[1]->contains(edge) || 
            faces[2]->contains(edge) || 
            faces[3]->contains(edge);
    }

    bool contains(const vertex_type* vert) const
    {
        return 
            faces[0]->contains(vert) ||
            faces[1]->contains(vert) ||
            faces[2]->contains(vert) ||
            faces[3]->contains(vert);
    }

    simplex(const simplex& simp)
    {
        faces = simp.faces;
    }
    simplex(const std::array<face_type*, 4>& faces)
        : faces(faces) {}
    simplex(const face_type* f0, const face_type* f1, const face_type* f2, const face_type* f3)
        : faces({ f0, f1, f2, f3 }) {}
};

template <std::size_t Dim, typename Real>
struct simplex<2, Dim, Real>
{
    static constexpr auto n = 2;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using facet_type = spt::edge<Dim, Real>;
    using edge_type = facet_type;
    using vertex_type = typename edge_type::vertex_type;

    std::array<edge_type*, 3> edges = { nullptr, };

    bool empty() const
    {
        return !edges[0] && !edges[1] && !edges[2];
    }

    template <typename SubPolytope>
    auto all_of() const
    {
        static_assert(std::is_same<edge_type, SubPolytope>() || std::is_same<vertex_type, SubPolytope>());
        if constexpr (std::is_same<edge_type, SubPolytope>())
            return edges;
        else
        {
            std::array<SubPolytope*, 4> verts = { nullptr, };
            std::size_t idx = 0;
            for (const auto& edge : edges)
                for (const auto& vert : edge->vertices)
                    if (std::find(verts.begin(), verts.begin() + idx, vert) == verts.begin() + idx)
                        verts[idx++] = vert;

            return verts;
        }
    }

    real_type area() const
    {
        vec<3> vec0 = edges[0]->vertices[1]->pos - edges[0]->vertices[0]->pos;
        vec<3> vec1 = edges[1]->vertices[1]->pos - edges[1]->vertices[0]->pos;
        return static_cast<real_type>(0.5) * vec<3>::cross(vec0, vec1).magnitude();
    }
    
    bool contains(const edge_type* edge) const
    {
        return edges[0] == edge || edges[1] == edge || edges[2] == edge;
    }

    bool contains(const vertex_type* vert) const
    {
        return edges[0]->contains(vert) || edges[1]->contains(vert) || edges[2]->contains(vert);
    }

    simplex(const simplex& simp)
    {
        edges = simp.edges;
    }
    simplex(const std::array<edge_type*, 3>& edges)
        : edges(edges) {}
    simplex(const edge_type* e0, const edge_type* e1, const edge_type* e2)
        : edges({ e0, e1, e2 }) {}
};

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
struct simplex_v
{
    static constexpr auto n = N;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, Real>;

    std::array<vertex_type*, N + 1> vertices;

    bool empty() const
    {
        return std::none_of(vertices.begin(), vertices.end(), [](auto x) -> bool { return x; });
    }

    bool contains(const vertex_type* vert) const
    {
        return std::find(vertices.begin(), vertices.end(), vert) != vertices.end();
    }

    simplex_v(const simplex_v& other)
    {
        vertices = other.vertices;
    }
    simplex_v(const std::array<vertex_type*, N + 1>& verts)
        : vertices(vertices) {}
    template <typename... Vertices>
    simplex_v(Vertices... verts)
        : vertices({ verts... }) {}
};

} // namespace spt
