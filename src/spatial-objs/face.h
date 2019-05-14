// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include "spatial-objs/edge.h"
#include "spatial-objs/vert.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Face
{
public:
    std::array<pmg::Edge*, 3> edges;

    vec3   computeCenter()  const;
    real_t computeQuality() const;
    real_t computeArea()    const;

//    static pmg::Edge* adjByEdge( const pmg::Face* face0, const pmg::Face* face1 );

    pmg::Vert* findVertNot( const pmg::Edge* edge ) const;
    pmg::Edge* findEdgeNot( const pmg::Vert* vert ) const;
    pmg::Edge* findEdge( const pmg::Vert* vert0, const pmg::Vert* vert1 ) const;
    pmg::Edge* shortestEdge() const;
    pmg::Edge* longestEdge()  const;

    bool contains( const pmg::Edge* edge ) const;
    bool contains( const pmg::Vert* vert ) const;

    Face( const pmg::Edge* edge0, const pmg::Edge* edge1, const pmg::Edge* edge2 );
    Face( const pmg::Vert* vert0, const pmg::Vert* vert1, const pmg::Vert* vert2 );
};

} // namespace pmg
