// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/edge.h"
#include "spatial-objs/vert.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Face
{
public:
    Edge* edges[3];

    Vec    computeCenter()  const;
    real_t computeQuality() const;
    real_t computeArea()    const;

    static Edge* intersectAlongEdge( const Face* facet0, const Face* facet1 );

    bool intersectsBy( const Vec& origin, const Vec& dir ) const;
    Vert* findVertNot( const Edge* edge ) const;
    Edge* findEdgeNot( const Vert* vert ) const;
    Edge* findEdge( const Vert* vert0, const Vert* vert1 ) const;
    Edge* findShortestEdge() const;
    Edge* findLongestEdge()  const;

    bool contains( const Edge*   edge ) const;
    bool contains( const Vert* vert ) const;

    Face( const Edge* edge0, const Edge* edge1, const Edge* edge2 );
    Face( const Vert* vert0, const Vert* vert1, const Vert* vert2 );
};

} // namespace pmg
