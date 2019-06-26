#pragma once
#include "mesh-base.h"
#include "../simplex.h"

namespace spt {

template <
    template <typename... Args> typename Pointer, 
    std::size_t N, std::size_t Dim, typename Real>
struct mesh_base<Pointer, spt::simplex<N, Dim, Real>> {
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, real_type>;
    // ...
    using facet_type = spt::simplex<N - 1, Dim, real_type>;
    using elem_type = spt::simplex<N, Dim, real_type>;

    std::vector<Pointer<vertex_type>> vertices;
    // ...
    std::vector<Pointer<facet_type>> facets;
    std::vector<Pointer<elem_type>> elements;

    mesh_base& operator=(const mesh_base& other) {
        vertices = other.vertices;
        // ...
        facets = other.facets;
        elements = other.elements;
    }
    mesh_base& operator=(mesh_base&& other) noexcept {
        vertices = std::move(other.vertices);
        // ...
        facets = std::move(other.facets);
        elements = std::move(other.elements);
    }

    mesh_base() {}
    mesh_base(const mesh_base& other) {
        *this = other;
    }
    mesh_base(mesh_base&& other) noexcept {
        *this = std::move(other);
    }
};


template <
    template <typename... Args> typename Pointer, 
    std::size_t N, std::size_t Dim, typename Real>
struct mesh_base<Pointer, spt::simplex_v<N, Dim, Real>> {
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, real_type>;
    using elem_type = spt::simplex_v<N, Dim, real_type>;

    std::vector<Pointer<vertex_type>> vertices;
    std::vector<Pointer<elem_type>> elements;

    mesh_base& operator=(const mesh_base& other) {
        vertices = other.vertices;
        elements = other.elements;
    }
    mesh_base& operator=(mesh_base&& other) noexcept {
        vertices = std::move(other.vertices);
        elements = std::move(other.elements);
    }

    mesh_base() {}
    mesh_base(const mesh_base& other) {
        *this = other;
    }
    mesh_base(mesh_base&& other) noexcept {
        *this = std::move(other);
    }
};

} // namespace pmg
