#pragma once
#include "dion.h"

namespace spt {

template <std::size_t Dim = 3, typename ValueType = typename spt::vec<Dim>::value_type>
using edge = dion<Dim, ValueType>;

} // namespace spt
