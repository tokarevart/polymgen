#pragma once
#include "../polyspt/simplex.h"


namespace pmg {

enum class elem_shape
{
    simplex,
    polytope
};

template <typename Polytope, elem_shape ElemShape, std::size_t N = Polytope::n>
struct indexed_mesh;

} // namespace pmg
