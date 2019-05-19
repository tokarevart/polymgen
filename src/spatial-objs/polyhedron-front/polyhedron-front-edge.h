// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/polyhedron.h"
#include "spatial-objs/polyhedron-front/polyhedron-front-face.h"
#include "spatial-objs/vert.h"
#include "spatial-objs/relations.h"
#include "real-type.h"

#include "definitions.h"

// TODO: maybe rename polyhedron-front-edge to volume-front-edge and pmg::front to pmg::volume::front
namespace pmg {
namespace front {

class Edge
{
    using pair_vv = std::pair<pmg::Vert*, pmg::Vert*>;
    using pair_ff = std::pair<front::Face*, front::Face*>;

public:
    pmg::Edge* edge;

    void refreshAngleData();

    real_t angle();
    real_t complexity();
    real_t computeAngle();
    real_t computeComplexity();

    // Opp means opposite.
    // NOTE: maybe store oppVert or information that there is no oppVert
    pair_vv      oppVerts();
    front::Edge* findOppEdge();

    // Adj means adjacent.
    pair_ff adjFFaces();

    bool addAdjFFace(       const front::Face* fFace );
    bool removeAdjFFace(    const front::Face* fFace );
    bool adjFFacesContains( const front::Face* fFace ) const;
    void fillAdjFFaces( const front::Face* fFace0, const front::Face* fFace1 );

    Edge( const Polyhedron* relatedPolyhedron, const pmg::Edge* edge );


private:
    Polyhedron* m_relatedPolyhedron;
    pair_ff m_adjFFaces = { nullptr, nullptr }; // TODO: make std::array

    real_t m_angle;
    real_t m_complexity;

    bool m_needAngleProcessing      = true;
    bool m_needComplexityProcessing = true;

    bool isAdjFacesFull();
    pair_ff fillAdjFFaces();
};

} // namespace front
} // namespace pmg
