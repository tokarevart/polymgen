#pragma once
#include "../polytope.h"
#include "genparams.h"
#include "../mesh/mesh.h"


namespace pmg {

template <typename Polytope, pmg::elem_shape ElemShape, std::size_t N = Polytope::n>
class mesher;

} // namespace pmg