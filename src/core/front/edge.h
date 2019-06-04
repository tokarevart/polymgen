// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "core/polyhedron.h"
#include "core/front/face.h"
#include "core/vert.h"
#include "core/relations.h"
#include "real-type.h"

#include "definitions.h"

// TODO: maybe rename front-edge to volume-front-edge and pmg::front to pmg::volume::front
namespace pmg::front {

class Edge
{
    using pair_vv = std::pair<pmg::Vert*, pmg::Vert*>;
    using pair_ff = std::pair<front::Face*, front::Face*>;

public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    pmg::Edge* edge;

    void refreshAngleData();

    real_t angle();
    real_t complexity();
    real_t computeAngle();
    real_t computeComplexity();

    // Opp means opposite.
    // TODO: rewrite these methods based on relations::opposite function
    // NOTE: maybe store oppVerts
    pair_vv      oppVerts();
    // NOTE: maybe store oppEdge and information of its existance
    front::Edge* findOppEdge();

    // Adj means adjacent.
    pair_ff adjFFaces();

    bool addAdjFFace(       const front::Face* fFace );
    bool removeAdjFFace(    const front::Face* fFace );
    bool adjFFacesContains( const front::Face* fFace ) const;
    void fillAdjFFaces( const front::Face* fFace0, const front::Face* fFace1 );
    // TODO: add method std::size_t nAdjFFaces() const;

    Edge( const Polyhedron* relatedPolyhedron, const pmg::Edge* edge );


private:
    Polyhedron* m_relatedPolyhedron; // TODO: it's not good to store it
    pair_ff m_adjFFaces = { nullptr, nullptr }; // TODO: make std::array instead

    real_t m_angle;
    real_t m_complexity;

    bool m_needAngleProcessing      = true;
    bool m_needComplexityProcessing = true;

    bool adjFacesFull();
    // TODO: remove this method and then before calling adjFFaces() adjacent faces must be assigned manually
    pair_ff fillAdjFFaces();
};

} // namespace pmg::front
