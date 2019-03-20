// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/polyhedron.h"
#include "spatial-objs/front/surface/front-surface-face.h"
#include "spatial-objs/vertex.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace front {
namespace surface {

class Edge
{
    using FrSuFace = front::surface::Face;
    using FrSuEdge  = front::surface::Edge;
    using pair_vv = std::pair<pmg::Vertex*, pmg::Vertex*>;
    using pair_ff = std::pair<FrSuFace*, FrSuFace*>;

public:
    pmg::Edge* edge;

    void refreshAngleData();

    real_t complexity();
    real_t angleExCos();
    real_t computeComplexity();
    real_t computeAngleExCos();
    real_t computeAngle();

    // Opp means opposite.
    pair_vv   findOppVerts();
    FrSuEdge* findOppEdge();

    // Adj means adjacent.
    pair_ff getAdjFFaces();

    bool addAdjFFace(       const FrSuFace* fFace );
    bool removeAdjFFace(    const FrSuFace* fFace );
    bool adjFFacesContains( const FrSuFace* fFace ) const;
    void fillAdjFFaces( const FrSuFace* fFace0, const FrSuFace* fFace1 );

    Edge( const Polyhedron* relatedPolyhedron, const pmg::Edge* edge );


private:
    Polyhedron* m_relatedPolyhedron;
    pair_ff m_adjFFaces{ nullptr, nullptr };

    real_t m_exCos;
    real_t m_complexity;

    bool m_needExCosProcessing      = true;
    bool m_needComplexityProcessing = true;

    bool isAdjFacesFull();
    pair_ff fillAdjFFaces();
};

} // namespace surface
} // namespace front
} // namespace pmg
