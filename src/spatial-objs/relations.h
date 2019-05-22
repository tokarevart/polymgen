#pragma once
#include "spatial-objs/polyhedron-front/polyhedron-front-face.h"
#include "spatial-objs/polyhedron-front/polyhedron-front-edge.h"
#include "spatial-objs/face.h"
#include "spatial-objs/edge.h"
#include "spatial-objs/vert.h"
#include "helpers/spatial-algs/vec.h"

namespace pmg {
namespace relations {

// TODO: replace with template <typename ByT> adjacent( ... ); and then template<> adjacent<each ByT> adjacent( ... );
pmg::Edge* adjacentByEdge( const pmg::Face* face0, const pmg::Face* face1 );
pmg::Vert* adjacentByVert( const pmg::Face* face0, const pmg::Face* face1 );
pmg::Vert* adjacentByVert( const pmg::Face* face,  const pmg::Edge* edge  );
pmg::Vert* adjacentByVert( const pmg::Edge* edge,  const pmg::Face* face  );
pmg::Vert* adjacentByVert( const pmg::Edge* edge0, const pmg::Edge* edge1 );

} // namespace relations
} // namespace pmg
