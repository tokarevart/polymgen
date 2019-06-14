#pragma once
#include "simplex.h"

namespace pmg {

enum class elem_shape
{
    simplex,
    polytope
};


template <typename Polytope, elem_shape ElemShape, std::size_t N = Polytope::n>
struct mesh;

template <typename Polytope, std::size_t N>
struct mesh<Polytope, elem_shape::simplex, N>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using vertex_type = spt::vertex<polytope_type::dim, real_type>;
    // ...
    using facet_type = spt::simplex<N - 1, polytope_type::dim, real_type>;
    using elem_type = spt::simplex<N, polytope_type::dim, real_type>;

    std::vector<vertex_type*> vertices;
    // ...
    std::vector<facet_type*> facets;
    std::vector<elem_type*> elements;
};

template <typename Polytope, std::size_t N>
struct mesh<Polytope, elem_shape::polytope, N>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using vertex_type = spt::vertex<polytope_type::dim, real_type>;
    // ...
    using facet_type = spt::polytope<N - 1, polytope_type::dim, real_type>;
    using elem_type = spt::polytope<N, polytope_type::dim, real_type>;

    std::vector<vertex_type*> vertices;
    // ...
    std::vector<facet_type*> facets;
    std::vector<elem_type*> elements;
};


template <typename Polytope, elem_shape ElemShape, std::size_t N>
struct mesh<spt::aggregate<Polytope*>, ElemShape, N>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using submesh_type = mesh<polytope_type, ElemShape, N>;

    std::vector<submesh_type*> meshes;
};

} // namespace pmg
