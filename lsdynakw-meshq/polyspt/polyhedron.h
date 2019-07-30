#pragma once
#include "face.h"

namespace spt {

template <std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::value_type>
using polyhedron = polytope<3, Dim, Real>;

template <std::size_t Dim, typename Real>
struct polytope<3, Dim, Real> {
    static constexpr std::size_t n = 3;
    static constexpr std::size_t dim = Dim;
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, Real>;
    using edge_type = spt::edge<Dim, Real>;
    using face_type = spt::face<Dim, Real>;
    using facet_type = face_type;

    std::vector<face_type*> faces;

    template <typename SubPolytope>
    auto all_of() const {
        static_assert(
            std::is_same<face_type, SubPolytope>() ||
            std::is_same<edge_type, SubPolytope>() ||
            std::is_same<vertex_type, SubPolytope>());
        if constexpr (std::is_same<face_type, SubPolytope>())
            return faces;
        else if constexpr (std::is_same<edge_type, SubPolytope>()) {
            std::vector<edge_type*> edges;
            for (const auto& face : faces)
                for (const auto& edge : face->edges)
                    if (std::find(edges.begin(), edges.begin() + edges.size(), edge) == edges.begin() + edges.size())
                        edges.push_back(edge);

            return edges;
        } else {
            std::vector<vertex_type*> verts;
            for (const auto& face : faces)
                for (const auto& edge : face->edges)
                    for (const auto& vert : edge->vertices)
                        if (std::find(verts.begin(), verts.begin() + verts.size(), vert) == verts.begin() + verts.size())
                            verts.push_back(vert);

            return verts;
        }
    }

    bool empty() const {
        return faces.empty();
    }
    
    bool contains(const face_type* face) const {
        return faces[0] == face || faces[1] == face || faces[2] == face;
    }
    bool contains(const edge_type* edge) const {
        return faces[0]->contains(edge) || faces[1]->contains(edge) || faces[2]->contains(edge);
    }
    bool contains(const vertex_type* vert) const {
        return faces[0]->contains(vert) || faces[1]->contains(vert) || faces[2]->contains(vert);
    }

    polytope(const polytope& poly) {
        faces = poly.faces;
    }
    polytope(polytope&& poly) {
        faces = std::move(poly.faces);
    }
    polytope(const std::vector<face_type*>& faces)
        : faces{ faces } {}
    template <typename... Faces>
    polytope(Faces... faces)
        : faces{ const_cast<face_type*>(faces)... } {}
};

} // namespace spt
