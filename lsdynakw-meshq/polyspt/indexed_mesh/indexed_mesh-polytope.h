#pragma once
#include <vector>
#include "indexed_mesh-base.h"
#include "../cell.h"


namespace spt {

template <std::size_t Dim, typename Real>
struct indexed_mesh<spt::edge<Dim, Real>> {
    using real_type = Real;
    using index_type = std::size_t;
    using vertex_type = spt::vec<Dim, real_type>;
    using edge_type = std::array<index_type, 2>;

    std::vector<vertex_type> vertices;
    std::vector<edge_type> edges;

    indexed_mesh() {}
    indexed_mesh(const indexed_mesh& other) {
        vertices = other.vertices;
        edges = other.edges;
    }
    indexed_mesh(indexed_mesh&& other) noexcept {
        vertices = std::move(other.vertices);
        edges = std::move(other.edges);
    }
};


template <std::size_t Dim, typename Real>
struct indexed_mesh<spt::face<Dim, Real>> {
    using real_type = Real;
    using index_type = std::size_t;
    using vertex_type = spt::vec<Dim, real_type>;
    using edge_type = std::array<index_type, 2>;
    using face_type = std::vector<index_type>;

    std::vector<vertex_type> vertices;
    std::vector<edge_type> edges;
    std::vector<face_type> faces;

    indexed_mesh() {}
    indexed_mesh(const indexed_mesh& other) {
        vertices = other.vertices;
        edges = other.edges;
        facets = other.facets;
    }
    indexed_mesh(indexed_mesh&& other) noexcept {
        vertices = std::move(other.vertices);
        edges = std::move(other.edges);
        facets = std::move(other.facets);
    }
};


template <std::size_t Dim, typename Real>
struct indexed_mesh<spt::cell<Dim, Real>> {
    using real_type = Real;
    using index_type = std::size_t;
    using vertex_type = spt::vec<Dim, real_type>;
    using edge_type = std::array<index_type, 2>;
    using face_type = std::vector<index_type>;
    using cell_type = std::vector<index_type>;

    std::vector<vertex_type> vertices;
    std::vector<edge_type> edges;
    std::vector<face_type> faces;
    std::vector<cell_type> cells;

    indexed_mesh() {}
    indexed_mesh(const indexed_mesh& other) {
        vertices = other.vertices;
        edges = other.edges;
        facets = other.facets;
        cells = other.cells;
    }
    indexed_mesh(indexed_mesh&& other) noexcept {
        vertices = std::move(other.vertices);
        edges = std::move(other.edges);
        facets = std::move(other.facets);
        cells = std::move(other.cells);
    }
};

} // namespace spt
