#pragma once
#include "mesh-base.h"

namespace spt {
       
template <typename Polytope>
struct mesh<Polytope, elem_shape::simplex>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using vertex_type = spt::vertex<polytope_type::dim, real_type>;
    // ...
    using facet_type = spt::simplex<polytope_type::n - 1, polytope_type::dim, real_type>;
    using elem_type = spt::simplex<polytope_type::n, polytope_type::dim, real_type>;

    std::vector<vertex_type*> vertices;
    // ...
    std::vector<facet_type*> facets;
    std::vector<elem_type*> elements;
};

template <typename Polytope>
struct mesh_v<Polytope, elem_shape::simplex>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using vertex_type = spt::vertex<polytope_type::dim, real_type>;
    using elem_type = spt::simplex_v<polytope_type::n, polytope_type::dim, real_type>;

    std::vector<vertex_type*> vertices;
    std::vector<elem_type*> elements;
};

} // namespace pmg
