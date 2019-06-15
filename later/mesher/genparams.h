#pragma once
#include "polytope.h"


namespace pmg {

template <typename Polytope>
struct genparams
{
    using polytope_type = Polytope;

    auto facet_genparams = genparams<typename polytope_type::facet_type>();

    // mesh generation parameters...
};

} // namespace pmg
