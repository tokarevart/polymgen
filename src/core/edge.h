// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include <list>
#include <vector>
#include "face.h"
#include "vert.h"
#include "../real-type.h"

#include "../definitions.h"


namespace pmg {

class Edge
{
    using pair_ff = std::pair<pmg::Face*, pmg::Face*>;

public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    std::array<pmg::Vert*, 2> verts;

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    pmg::Vert* findNot( const pmg::Edge* edge ) const;
    pmg::Vert* findNot( const pmg::Vert* vert ) const;

    // TODO: move these methods to some file like surface/relations.h or something like that
    void flip(         std::list<pmg::Edge*>& edgesList, std::list<pmg::Face*>& facesList );
    bool flipIfNeeded( std::list<pmg::Edge*>& edgesList, std::list<pmg::Face*>& facesList );
    // Adj means adjacent.
    // TODO: move these methods to pmg::relations
    std::list<pmg::Face*> findAdjFaces(  const std::list<pmg::Face*>& facesList ) const;
    pair_ff               find2AdjFaces( const std::list<pmg::Face*>& facesList ) const;

    bool contains( const pmg::Vert* vert ) const;

    // TODO: move this method to pmg::relations
    bool belongsToShell();

    bool needToFlip( const std::list<pmg::Face*>& facesList );

    Edge( const pmg::Vert* vert0, const pmg::Vert* vert1 );
};

} // namespace pmg
