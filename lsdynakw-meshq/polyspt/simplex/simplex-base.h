#pragma once
#include <cstddef>
#include <array>
#include <algorithm>
#include "../face.h"

namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename ValueType = typename spt::vec<Dim>::value_type>
struct simplex;

template <std::size_t N, std::size_t Dim = 3, typename ValueType = typename spt::vec<Dim>::value_type>
struct simplex_v;

} // namespace spt
