// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include "edge.h"
#include "vert.h"
#include "../helpers/spatial/vec.h"
#include "../real-type.h"

#include "../definitions.h"


namespace pmg {

// TODO: generalize geom structures with templates (and no interfaces); use Polytope and so on; with static and dynamic num of subpolytopes
class Face {
public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    std::array<pmg::Edge*, 3> edges;

    // TODO: add method std::array<pmg::Vert*, 3> verts() const;

    spt::vec3 computeCenter() const;
    real_t computeQuality() const;
    real_t computeArea() const;

//    static pmg::Edge* adjByEdge( const pmg::Face* face0, const pmg::Face* face1 );

    pmg::Vert* findVertNot(const pmg::Edge* edge) const;
    pmg::Edge* findEdgeNot(const pmg::Vert* vert) const;
    pmg::Edge* findEdge(const pmg::Vert* vert0, const pmg::Vert* vert1) const;
    pmg::Edge* shortestEdge() const;
    pmg::Edge* longestEdge() const;

    bool contains(const pmg::Edge* edge) const;
    bool contains(const pmg::Vert* vert) const;

    Face(const pmg::Edge* edge0, const pmg::Edge* edge1, const pmg::Edge* edge2);
    Face(const pmg::Vert* vert0, const pmg::Vert* vert1, const pmg::Vert* vert2);
};

} // namespace pmg
