#pragma once
#include "genparams.h"
#include "../polyspt/mesh.h"


namespace pmg {

template <typename Polytope, spt::elem_shape ElemShape, std::size_t N = Polytope::n>
class mesher;

} // namespace pmg