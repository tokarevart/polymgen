#pragma once
#include <vector>
#include <memory>
#include "../polytope-base.h"

namespace spt {

template <typename T>
using raw_ptr = T*;

template <
    template <typename ElemType, typename... Args> typename Pointer, 
    typename ElemType, 
    template <typename Polytope> typename HowMuchPolytopes = spt::single>
struct mesh_base;

template <typename ElemType, template <typename Polytope> typename HowMuchPolytopes = spt::single>
using unique_mesh = mesh_base<std::unique_ptr, ElemType, HowMuchPolytopes>;

template <typename ElemType, template <typename Polytope> typename HowMuchPolytopes = spt::single>
using raw_mesh = mesh_base<raw_ptr, ElemType, HowMuchPolytopes>;

template <template <typename MeshType, typename... Args> typename Pointer, typename ElemType>
struct mesh_base<Pointer, ElemType, spt::composition> {
    std::vector<Pointer<mesh_base<Pointer, ElemType>>> meshes;
};

} // namespace pmg
