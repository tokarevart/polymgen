#pragma once
#include "mesh-base.h"

namespace pmg {
       
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

} // namespace pmg
