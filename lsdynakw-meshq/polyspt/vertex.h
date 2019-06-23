#pragma once
#include "monon.h"

namespace spt {

template <std::size_t Dim = 3, typename ValueType = typename spt::vec<Dim>::value_type>
using vertex = monon<Dim, ValueType>;

} // namespace spt
