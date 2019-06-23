#pragma once
#include <cstddef>
#include "vec.h"


namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
struct polytope;


template <typename Polytope>
struct aggregate;

} // namespace spt
