#pragma once
#include "mesh-base.h"

namespace spt {

template <
    typename ElemType,
    spt::amount MeshesAmount = spt::amount::single>
struct make_unique_mesh_t;

template <
    typename ElemType,
    spt::amount MeshesAmount = spt::amount::single>
struct make_raw_mesh_t;


template <typename ElemType>
struct make_unique_mesh_t<ElemType, spt::amount::single> {
    template <template <typename... Args> typename Pointer> 
    static spt::unique_mesh<ElemType> 
        run(const spt::mesh_base<Pointer, ElemType>& mesh) {

    }
};

template <typename ElemType>
struct make_raw_mesh_t<ElemType, spt::amount::single> {
    template <template <typename... Args> typename Pointer>
    static spt::raw_mesh<ElemType> 
        run(const spt::mesh_base<Pointer, ElemType>& mesh) {

    }
};

template <typename ElemType>
struct make_unique_mesh_t<ElemType, spt::amount::aggregate> {
    template <template <typename... Args> typename Pointer>
    static spt::unique_mesh<ElemType, spt::amount::aggregate>
        run(const spt::mesh_base<Pointer, ElemType, spt::amount::aggregate>& mesh) {

    }
};

template <typename ElemType>
struct make_raw_mesh_t<ElemType, spt::amount::aggregate> {
    template <template <typename... Args> typename Pointer>
    static spt::raw_mesh<ElemType, spt::amount::aggregate>
        run(const spt::mesh_base<Pointer, ElemType, spt::amount::aggregate>& mesh) {

    }
};

template <
    typename ElemType, 
    spt::amount MeshesAmount = spt::amount::single, 
    template <typename... Args> typename Pointer>
auto make_unique_mesh(const spt::mesh_base<Pointer, ElemType, MeshesAmount>& mesh) {
    return make_unique_mesh_t<ElemType, MeshesAmount>::run(mesh);
}

template <
    typename ElemType, 
    spt::amount MeshesAmount = spt::amount::single, 
    template <typename... Args> typename Pointer>
auto make_raw_mesh(const spt::mesh_base<Pointer, ElemType, MeshesAmount>& mesh) {
    return make_raw_mesh_t<ElemType, MeshesAmount>::run(mesh);
}

} // namespace spt