#pragma once


namespace pmg {

enum class elem_shape
{
    simplex,
    polytope
};

template <typename Polytope, elem_shape ElemShape>
struct indexed_mesh;

} // namespace pmg
