#pragma once
#include "../polyspt/simplex.h"

namespace pmg {

enum class elem_shape
{
    simplex,
    polytope
};


template <typename Polytope, elem_shape ElemShape>
struct mesh;

template <typename Polytope, elem_shape ElemShape>
struct mesh<spt::aggregate<Polytope>, ElemShape>
{
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;
    using submesh_type = mesh<polytope_type, ElemShape>;

    std::vector<submesh_type*> meshes;
};

} // namespace pmg
