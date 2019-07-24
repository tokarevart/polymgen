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
    spt::amount MeshesAmount = spt::amount::single>
struct mesh_base; // TODO: maybe use tuple inheritance or composition or something like that

template <typename ElemType, spt::amount MeshesAmount = spt::amount::single>
using unique_mesh = mesh_base<std::unique_ptr, ElemType, MeshesAmount>;

template <typename ElemType, spt::amount MeshesAmount = spt::amount::single>
using raw_mesh = mesh_base<raw_ptr, ElemType, MeshesAmount>;

template <template <typename... Args> typename Pointer, typename ElemType>
struct mesh_base<Pointer, ElemType, spt::amount::aggregate> {
    std::vector<Pointer<mesh_base<Pointer, ElemType>>> meshes;

    template <typename = std::enable_if_t<std::is_same_v<std::unique_ptr<void>, Pointer<void>>>>
    auto get() { return spt::to_raw_mesh(*this); }
};

} // namespace spt
