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
    struct to_raw_mesh_t;


template <typename ElemType>
struct make_unique_mesh_t<ElemType, spt::amount::single> {
    static spt::unique_mesh<ElemType> run(const spt::raw_mesh<ElemType>& mesh) {

    }
};

template <typename ElemType>
struct to_raw_mesh_t<ElemType, spt::amount::single> {
    static spt::raw_mesh<ElemType> run(const spt::unique_mesh<ElemType>& mesh) {

    }
};

template <typename ElemType>
struct make_unique_mesh_t<ElemType, spt::amount::aggregate> {
    static spt::unique_mesh<ElemType, spt::amount::aggregate>
        run(const spt::raw_mesh<ElemType, spt::amount::aggregate>& mesh) {

    }
};

template <typename ElemType>
struct to_raw_mesh_t<ElemType, spt::amount::aggregate> {
    static spt::raw_mesh<ElemType, spt::amount::aggregate>
        run(const spt::unique_mesh<ElemType, spt::amount::aggregate>& mesh) {

    }
};

template <typename ElemType, spt::amount MeshesAmount = spt::amount::single>
auto make_unique_mesh(const spt::raw_mesh<ElemType, MeshesAmount>& mesh) {
    return make_unique_mesh_t<ElemType, MeshesAmount>::run(mesh);
}

template <typename ElemType, spt::amount MeshesAmount = spt::amount::single>
auto to_raw_mesh(const spt::unique_mesh<ElemType, MeshesAmount>& mesh) {
    return to_raw_mesh_t<ElemType, MeshesAmount>::run(mesh);
}

} // namespace spt