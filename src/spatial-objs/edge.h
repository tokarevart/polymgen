#pragma once
#include <list>
#include <vector>
#include "spatial-objs/facet.h"
#include "spatial-objs/vertex.h"

#include "definitions.h"


namespace pmg {

class Edge
{
    using pair_ff = std::pair<pmg::Facet*, pmg::Facet*>;

public:
    Vertex* verts[2];

    double magnitude()    const;
    double sqrMagnitude() const;

    Vertex* findNot( const Edge*   edge ) const;
    Vertex* findNot( const Vertex* vert ) const;

    void flip( std::list<Edge*>& edgesList, std::list<Facet*>& facetsList );
    bool flipIfNeeded( std::list<Edge*>& edgesList, std::list<Facet*>& facetsList );
    // Adj means adjacent.
    void findAdjFacets( const std::list<Facet*>& facetsList, std::list<Facet*>& adjFacets ) const;
    pair_ff find2AdjFacets( const std::list<Facet*>& facetsList ) const;

    bool contains( const Vertex* vert ) const;

    bool belongsToShell();

    bool needToFlip( const std::list<Facet*>& facetsList );

    Edge( const Vertex* vert0, const Vertex* vert1 );
};

} // namespace pmg
