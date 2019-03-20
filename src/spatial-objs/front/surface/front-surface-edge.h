// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/crystallite.h"
#include "spatial-objs/front/surface/front-surface-facet.h"
#include "spatial-objs/vertex.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace front {
namespace surface {

class Edge
{
    using FrSuFacet = front::surface::Facet;
    using FrSuEdge  = front::surface::Edge;
    using pair_vv = std::pair<pmg::Vertex*, pmg::Vertex*>;
    using pair_ff = std::pair<FrSuFacet*, FrSuFacet*>;

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
    pair_ff getAdjFFacets();

    bool addAdjFFacet(       const FrSuFacet* fFacet );
    bool removeAdjFFacet(    const FrSuFacet* fFacet );
    bool adjFFacetsContains( const FrSuFacet* fFacet ) const;
    void fillAdjFFacets( const FrSuFacet* fFacet0, const FrSuFacet* fFacet1 );

    Edge( const Crystallite* relatedCrys, const pmg::Edge* edge );


private:
    Crystallite* m_relatedCrys;
    pair_ff m_adjFFacets{ nullptr, nullptr };

    real_t m_exCos;
    real_t m_complexity;

    bool m_needExCosProcessing      = true;
    bool m_needComplexityProcessing = true;

    bool isAdjFacetsFull();
    pair_ff fillAdjFFacets();
};

} // namespace surface
} // namespace front
} // namespace pmg
