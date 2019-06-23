#pragma once
#include <cstddef>
#include "vec.h"


namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename ValueType = typename spt::vec<Dim>::value_type>
struct polytope;


template <typename Polytope>
struct aggregate;

} // namespace spt
