#pragma once
#include "mesh-base.h"

namespace spt {

template <std::size_t N, std::size_t Dim, typename Real>
struct mesh_base<std::unique_ptr<spt::polytope<N, Dim, Real>>> {
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, real_type>;
    // ...
    using facet_type = spt::polytope<N - 1, Dim, real_type>;
    using elem_type = spt::polytope<N, Dim, real_type>;

    std::vector<std::unique_ptr<vertex_type>> vertices;
    // ...
    std::vector<std::unique_ptr<facet_type>> facets;
    std::vector<std::unique_ptr<elem_type>> elements;

    mesh_base() {}
    mesh_base(mesh_base&& other) noexcept {
        vertices = std::move(other.vertices);
        // ...
        facets = std::move(other.facets);
        elements = std::move(other.elements);
    }
};


template <std::size_t N, std::size_t Dim, typename Real>
struct mesh_base<spt::polytope<N, Dim, Real>*> {
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, real_type>;
    // ...
    using facet_type = spt::polytope<N - 1, Dim, real_type>;
    using elem_type = spt::polytope<N, Dim, real_type>;

    std::vector<vertex_type*> vertices;
    // ...
    std::vector<facet_type*> facets;
    std::vector<elem_type*> elements;

    mesh_base() {}
    mesh_base(const mesh_base& other) {
        vertices = other.vertices;
        // ...
        facets = other.facets;
        elements = other.elements;
    }
    mesh_base(mesh_base&& other) noexcept {
        vertices = std::move(other.vertices);
        // ...
        facets = std::move(other.facets);
        elements = std::move(other.elements);
    }
};

} // namespace pmg
