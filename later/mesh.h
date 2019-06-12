#pragma once
#include "polytope.h"

namespace pmg {

template <typename Polytope>
struct mesh
{
    using polytope_type = Polytope;

    // ...
};

template <typename Polytope>
struct mesh<std::vector<Polytope>>
{
    using polytope_type = Polytope;

    // ...
};

template <typename Polytope>
struct raw_mesh
{
    using polytope_type = Polytope;

    // ...
};

template <typename Polytope>
struct raw_mesh<std::vector<Polytope>>
{
    using polytope_type = Polytope;

    // ...
};

} // namespace pmg
