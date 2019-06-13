#pragma once
#include "simplex.h"

namespace pmg {

template <typename Polytope>
struct indexed_mesh
{
    using polytope_type = Polytope;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, typename polytope_type::real_type>;
    using simplex_type = std::array<index_type, polytope_type::n + 1>;

    std::vector<vertex_type> vertices;
    std::vector<simplex_type> simplices;
};


template <typename Polytope>
struct indexed_mesh<spt::aggregate<Polytope>>
{
    using polytope_type = Polytope;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, typename polytope_type::real_type>;
    using simplex_type = std::array<index_type, polytope_type::n + 1>;
    using submesh_type = std::vector<simplex_type>;
    
    std::vector<vertex_type> vertices;
    std::vector<submesh_type> meshes;
};


template <typename Polytope>
struct mesh
{
    using polytope_type = Polytope;
    using vertex_type = spt::vertex<polytope_type::dim, typename polytope_type::real_type>;
    // ...
    using facet_type = spt::simplex<polytope_type::n - 1, polytope_type::dim, typename polytope_type::real_type>;
    using simplex_type = spt::simplex<polytope_type::n, polytope_type::dim, typename polytope_type::real_type>;

    std::vector<vertex_type*> vertices;
    // ...
    std::vector<facet_type*> facets;
    std::vector<simplex_type*> simplices;
};


template <typename Polytope>
struct mesh<spt::aggregate<Polytope>>
{
    using polytope_type = Polytope;
    using submesh_type = mesh<polytope_type>;

    std::vector<submesh_type*> meshes;
};

} // namespace pmg
