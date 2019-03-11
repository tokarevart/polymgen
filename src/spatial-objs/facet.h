#pragma once
#include "spatial-objs/edge.h"
#include "spatial-objs/vertex.h"
#include "helpers/spatial-algs/vec.h"

#include "definitions.h"


namespace pmg {

class Facet
{
    typedef tva::Vec Vec;
    typedef tva::Point Point;

public:
    Edge* edges[3];

    Vec    computeCenter()  const;
    double computeQuality() const;
    double computeArea()    const;

    static Edge* intersectAlongEdge( const Facet* facet0, const Facet* facet1 );

    bool intersectsBy( const Point& origin, const Vec& dir )     const;
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
