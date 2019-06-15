#pragma once
#include "indexed_mesh-base.h"


namespace pmg {

template <typename Polytope>
struct indexed_mesh<Polytope, pmg::elem_shape::polytope, 1>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, real_type>;
    using elem_type = std::array<index_type, 2>;

    std::vector<vertex_type> vertices;
    std::vector<elem_type> elements;
};

template <typename Polytope>
struct indexed_mesh<Polytope, pmg::elem_shape::polytope, 2>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, real_type>;
    using edge_type = std::array<index_type, 2>;
    using elem_type = std::vector<index_type>;

    std::vector<vertex_type> vertices;
    std::vector<edge_type> edges;
    std::vector<elem_type> elements;
};

template <typename Polytope>
struct indexed_mesh<Polytope, pmg::elem_shape::polytope, 3>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using index_type = std::size_t;
    using vertex_type = spt::vec<polytope_type::dim, real_type>;
    using edge_type = std::array<index_type, 2>;
    using face_type = std::vector<index_type>;
    using elem_type = std::vector<index_type>;

    std::vector<vertex_type> vertices;
    std::vector<edge_type> edges;
    std::vector<face_type> faces;
    std::vector<elem_type> elements;
};

} // namespace pmg
