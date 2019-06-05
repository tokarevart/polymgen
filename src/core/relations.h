#pragma once
#include "front/face.h"
#include "front/edge.h"
#include "face.h"
#include "edge.h"
#include "vert.h"
#include "../helpers/spatial/vec.h"

namespace pmg {
namespace relations {

// TODO: replace with template <typename By> By* adjacent( ... ); and then template<> By* adjacent<each By> adjacent( ... );
pmg::Edge* adjacentByEdge( const pmg::Face* face0, const pmg::Face* face1 );
pmg::Vert* adjacentByVert( const pmg::Face* face0, const pmg::Face* face1 );
pmg::Vert* adjacentByVert( const pmg::Face* face,  const pmg::Edge* edge  );
pmg::Vert* adjacentByVert( const pmg::Edge* edge,  const pmg::Face* face  );
pmg::Vert* adjacentByVert( const pmg::Edge* edge0, const pmg::Edge* edge1 );
// TODO: place all the functions related to front, surface and volume in the corresponding files (NOT HERE)

template <typename Opp> std::array<Opp*, 2> opposite( const pmg::Edge* edge ); // or {}, i don't remember
template <typename Opp> Opp*                opposite( const pmg::Edge* edge );

template<> std::array<pmg::Vert*, 2> opposite( const pmg::Edge* edge );
template<> pmg::Edge*                opposite( const pmg::Edge* edge );

// TODO: add template partial specialization for adjacent( ... ) to find adjacent front faces to front edge

} // namespace relations
} // namespace pmg
