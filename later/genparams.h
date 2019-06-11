#pragma once
#include "polytope.h"


namespace pmg {

template <typename Polytope>
struct genparams
{
    using polytope_type = Polytope;

    auto facet_genparams = genparams<spt::polytope<N-1, Dim, Real>>();

    // mesh generation parameters...
};

} // namespace pmg
