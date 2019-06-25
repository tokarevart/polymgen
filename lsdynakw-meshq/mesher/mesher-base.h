#pragma once
#include "genparams.h"
#include "../polyspt/mesh.h"


namespace pmg {

template <
    typename Polytope, 
    template <std::size_t N, std::size_t Dim, typename Real> typename ElemType,
    template <typename... Args> typename Pointer>
class mesher;

} // namespace pmg