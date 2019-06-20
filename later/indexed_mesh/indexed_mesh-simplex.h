#pragma once
#include <vector>
#include "indexed_mesh-base.h"


namespace pmg {

template <typename Polytope>
struct indexed_mesh<Polytope, pmg::elem_shape::simplex>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, real_type>;
    using elem_type = std::array<index_type, polytope_type::n + 1>;

    std::vector<vertex_type> vertices;
    std::vector<elem_type> elements;
};

template <typename Polytope>
struct indexed_mesh<spt::aggregate<Polytope>, pmg::elem_shape::simplex>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, real_type>;
    using elem_type = std::array<index_type, polytope_type::n + 1>;
    using submesh_type = std::vector<elem_type>;

    std::vector<vertex_type> vertices;
    std::vector<submesh_type> meshes;
};

} // namespace pmg
