#pragma once
#include "spatial-objs/face.h"
#include "spatial-objs/edge.h"
#include "spatial-objs/vert.h"
#include "helpers/spatial-algs/vec.h"

namespace pmg {
namespace relations {

pmg::Edge* adjacentByEdge( const pmg::Face* face0, const pmg::Face* face1 );
pmg::Vert* adjacentByVert( const pmg::Face* face0, const pmg::Face* face1 );
pmg::Vert* adjacentByVert( const pmg::Face* face,  const pmg::Edge* edge  );
pmg::Vert* adjacentByVert( const pmg::Edge* edge,  const pmg::Face* face  );
pmg::Vert* adjacentByVert( const pmg::Edge* edge0, const pmg::Edge* edge1 );

bool contains( const front::Face* fFace,  const front::Edge* fEdge );
bool contains( const front::Face* fFace,  const   pmg::Edge* edge  );
bool contains( const front::Face* fFace,  const   pmg::Vert* vert  );
bool contains( const shell::Face* sFace,  const shell::Edge* sEdge );
bool contains( const shell::Face* sFace,  const shell::Vert* sVert );
bool contains( const shell::Face* sFace,  const   pmg::Face* face  );
bool contains( const shell::Face* sFace,  const   pmg::Edge* edge  );
bool contains( const shell::Face* sFace,  const   pmg::Vert* vert  );
bool contains( const shell::Edge* sEdge,  const shell::Vert* sVert );
bool contains( const shell::Edge* sEdge,  const   pmg::Edge* edge  );
bool contains( const shell::Edge* sEdge,  const   pmg::Vert* vert  );
bool contains( const   pmg::Face* face,   const   pmg::Edge* edge  );
bool contains( const   pmg::Face* face,   const   pmg::Vert* vert  );
bool contains( const   pmg::Edge* edge,   const   pmg::Vert* vert  );

} // namespace relations
} // namespace pmg
