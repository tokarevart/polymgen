// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/face.h"
#include "spatial-objs/vert.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Edge
{
    using pair_ff = std::pair<pmg::Face*, pmg::Face*>;

public:
    Vert* verts[2];

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    Vert* findNot( const Edge* edge ) const;
    Vert* findNot( const Vert* vert ) const;

    void flip( std::list<Edge*>& edgesList, std::list<Face*>& facesList );
    bool flipIfNeeded( std::list<Edge*>& edgesList, std::list<Face*>& facesList );
    // Adj means adjacent.
    void findAdjFaces( const std::list<Face*>& facesList, std::list<Face*>& adjFaces ) const;
    pair_ff find2AdjFaces( const std::list<Face*>& facesList ) const;

    bool contains( const Vert* vert ) const;

    bool belongsToShell();

    bool needToFlip( const std::list<Face*>& facesList );

    Edge( const Vert* vert0, const Vert* vert1 );
};

} // namespace pmg
