#pragma once
#include <vector>
#include "indexed_mesh-base.h"
#include "../simplex.h"


namespace spt {

template <std::size_t N, std::size_t Dim, typename Real>
struct indexed_mesh<spt::simplex<N, Dim, Real>> {
    using real_type = Real;
    using index_type = std::size_t;
    using vertex_type = spt::vec<Dim, real_type>;
    using elem_type = std::array<index_type, N + 1>;

    std::vector<vertex_type> vertices;
    std::vector<elem_type> elements;

    indexed_mesh() {}
    indexed_mesh(const indexed_mesh& other) {
        vertices = other.vertices;
        elements = other.elements;
    }
    indexed_mesh(indexed_mesh&& other) noexcept {
        vertices = std::move(other.vertices);
        elements = std::move(other.elements);
    }
};


template <std::size_t N, std::size_t Dim, typename Real>
struct indexed_mesh_composition<spt::simplex<N, Dim, Real>> {
    using real_type = Real;
    using index_type = std::size_t;
    using vertex_type = spt::vec<Dim, real_type>;
    using elem_type = std::array<index_type, N + 1>;
    using mesh_type = std::vector<elem_type>;

    std::vector<vertex_type> vertices;
    std::vector<mesh_type> meshes;

    indexed_mesh_composition() {}
    indexed_mesh_composition(const indexed_mesh_composition& other) {
        vertices = other.vertices;
        meshes = other.meshes;
    }
    indexed_mesh_composition(indexed_mesh_composition&& other) noexcept {
        vertices = std::move(other.vertices);
        meshes = std::move(other.meshes);
    }
};

} // namespace spt
