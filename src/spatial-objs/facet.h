// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/edge.h"
#include "spatial-objs/vertex.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Facet
{
public:
    Edge* edges[3];

    Vec    computeCenter()  const;
    real_t computeQuality() const;
    real_t computeArea()    const;

    static Edge* intersectAlongEdge( const Facet* facet0, const Facet* facet1 );

    bool intersectsBy( const Vec& origin, const Vec& dir )     const;
    Vertex* findVertNot( const Edge* edge )        const;
    Edge*   findEdgeNot( const Vertex* vert )      const;
    Edge*   findEdge( const Vertex* vert0, const Vertex* vert1 ) const;
    Edge* findShortestEdge() const;
    Edge* findLongestEdge()  const;

    bool contains( const Edge*   edge ) const;
    bool contains( const Vertex* vert ) const;

    Facet( const Edge*   edge0, const Edge*   edge1, const Edge*   edge2 );
    Facet( const Vertex* vert0, const Vertex* vert1, const Vertex* vert2 );
};

} // namespace pmg
