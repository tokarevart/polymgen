// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
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
    std::array<pmg::Vert*, 2> verts;

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    pmg::Vert* findNot( const pmg::Edge* edge ) const;
    pmg::Vert* findNot( const pmg::Vert* vert ) const;

    void flip(         std::list<Edge*>& edgesList, std::list<pmg::Face*>& facesList );
    bool flipIfNeeded( std::list<Edge*>& edgesList, std::list<pmg::Face*>& facesList );
    // Adj means adjacent.
    std::list<pmg::Face*> findAdjFaces(  const std::list<pmg::Face*>& facesList ) const;
    pair_ff               find2AdjFaces( const std::list<pmg::Face*>& facesList ) const;

    bool contains( const pmg::Vert* vert ) const;

    bool belongsToShell();

    bool needToFlip( const std::list<pmg::Face*>& facesList );

    Edge( const pmg::Vert* vert0, const pmg::Vert* vert1 );
};

} // namespace pmg
