#pragma once
#include <cstddef>
#include <array>
#include <list>
#include "vec.h"


namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
struct polytope;

template <typename Polytope>
using aggregate = std::list<Polytope*>;

} // namespace spt
