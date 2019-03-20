// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/face.h"
#include "spatial-objs/vertex.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Edge
{
    using pair_ff = std::pair<pmg::Face*, pmg::Face*>;

public:
    Vertex* verts[2];

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    Vertex* findNot( const Edge*   edge ) const;
    Vertex* findNot( const Vertex* vert ) const;

    void flip( std::list<Edge*>& edgesList, std::list<Face*>& facesList );
    bool flipIfNeeded( std::list<Edge*>& edgesList, std::list<Face*>& facesList );
    // Adj means adjacent.
    void findAdjFaces( const std::list<Face*>& facesList, std::list<Face*>& adjFaces ) const;
    pair_ff find2AdjFaces( const std::list<Face*>& facesList ) const;

    bool contains( const Vertex* vert ) const;

    bool belongsToShell();

    bool needToFlip( const std::list<Face*>& facesList );

    Edge( const Vertex* vert0, const Vertex* vert1 );
};

} // namespace pmg
