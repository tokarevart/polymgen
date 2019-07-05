#pragma once
#include <vector>
#include <memory>
#include "../polytope-base.h"

namespace spt {

template <typename T>
using raw_ptr = T*;

template <
    template <typename... Args> typename Pointer, 
    typename ElemType, 
    template <typename Polytope> typename HowMany = spt::single>
struct mesh_base; // TODO: maybe use tuple inheritance or composition or something like that

template <typename ElemType, template <typename Polytope> typename HowMany = spt::single>
using unique_mesh = mesh_base<std::unique_ptr, ElemType, HowMany>;

template <typename ElemType, template <typename Polytope> typename HowMany = spt::single>
using raw_mesh = mesh_base<raw_ptr, ElemType, HowMany>;

template <template <typename... Args> typename Pointer, typename ElemType>
struct mesh_base<Pointer, ElemType, spt::multi> {
    std::vector<Pointer<mesh_base<Pointer, ElemType>>> meshes;
};

} // namespace pmg
