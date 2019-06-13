#pragma once
#include "simplex.h"

namespace pmg {

enum class elem_shape
{
    simplex,
    polytope
};


template <typename Polytope, elem_shape ElemShape, std::size_t N = Polytope::n>
struct indexed_mesh;


template <typename Polytope, std::size_t N>
struct indexed_mesh<Polytope, elem_shape::simplex, N>
{
    using polytope_type = Polytope;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, typename polytope_type::real_type>;
    using elem_type = std::array<index_type, polytope_type::n + 1>;

    std::vector<vertex_type> vertices;
    std::vector<elem_type> elements;
};

template <typename Polytope>
struct indexed_mesh<Polytope, elem_shape::polytope, 1>
{
    using polytope_type = Polytope;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, typename polytope_type::real_type>;
    using elem_type = std::array<index_type, 2>;

    std::vector<vertex_type> vertices;
    std::vector<elem_type> elements;
};

template <typename Polytope>
struct indexed_mesh<Polytope, elem_shape::polytope, 2>
{
    using polytope_type = Polytope;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, typename polytope_type::real_type>;
    using edge_type = std::array<index_type, 2>;
    using elem_type = std::vector<index_type>;

    std::vector<vertex_type> vertices;
    std::vector<edge_type> edges;
    std::vector<elem_type> elements;
};

template <typename Polytope>
struct indexed_mesh<Polytope, elem_shape::polytope, 3>
{
    using polytope_type = Polytope;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, typename polytope_type::real_type>;
    using edge_type = std::array<index_type, 2>;
    using face_type = std::vector<index_type>;
    using elem_type = std::vector<index_type>;

    std::vector<vertex_type> vertices;
    std::vector<edge_type> edges;
    std::vector<face_type> faces;
    std::vector<elem_type> elements;
};


template <typename Polytope, std::size_t N>
struct indexed_mesh<spt::aggregate<Polytope>, elem_shape::simplex, N>
{
    using polytope_type = Polytope;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, typename polytope_type::real_type>;
    using simplex_type = std::array<index_type, polytope_type::n + 1>;
    using submesh_type = std::vector<simplex_type>;
    
    std::vector<vertex_type> vertices;
    std::vector<submesh_type> meshes;
};


template <typename Polytope, elem_shape ElemShape, std::size_t N = Polytope::n>
struct mesh;

template <typename Polytope, std::size_t N>
struct mesh<Polytope, elem_shape::simplex, N>
{
    using polytope_type = Polytope;
    using vertex_type = spt::vertex<polytope_type::dim, typename polytope_type::real_type>;
    // ...
    using facet_type = spt::simplex<polytope_type::n - 1, polytope_type::dim, typename polytope_type::real_type>;
    using elem_type = spt::simplex<polytope_type::n, polytope_type::dim, typename polytope_type::real_type>;

    std::vector<vertex_type*> vertices;
    // ...
    std::vector<facet_type*> facets;
    std::vector<elem_type*> elements;
};

template <typename Polytope, std::size_t N>
struct mesh<Polytope, elem_shape::polytope, N>
{
    using polytope_type = Polytope;
    using vertex_type = spt::vertex<polytope_type::dim, typename polytope_type::real_type>;
    // ...
    using facet_type = spt::polytope<polytope_type::n - 1, polytope_type::dim, typename polytope_type::real_type>;
    using elem_type = spt::polytope<polytope_type::n, polytope_type::dim, typename polytope_type::real_type>;

    std::vector<vertex_type*> vertices;
    // ...
    std::vector<facet_type*> facets;
    std::vector<elem_type*> elements;
};


template <typename Polytope, elem_shape ElemShape, std::size_t N>
struct mesh<spt::aggregate<Polytope>, ElemShape, N>
{
    using polytope_type = Polytope;
    using submesh_type = mesh<polytope_type, ElemShape>;

    std::vector<submesh_type*> meshes;
};

} // namespace pmg
