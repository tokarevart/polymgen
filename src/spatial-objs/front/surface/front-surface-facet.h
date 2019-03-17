// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/crystallite.h"
#include "helpers/spatial-algs/vec.h"

#include "definitions.h"


namespace pmg {
namespace front {
namespace surface {

class Facet
{
    using FrSuEdge = front::surface::Edge;
    using Vec = tva::Vec;

public:
    pmg::Facet* facet;
    FrSuEdge* fEdges[3] { nullptr, nullptr, nullptr };
    Vec normal;

    Vec    computeNormal();
    Vec    computeCenter();
    double computeQuality();

    FrSuEdge* findFEdge( const pmg::Edge* edge ) const;
    FrSuEdge* findFEdge( const pmg::Vertex* v0, const pmg::Vertex* v1 ) const;
    FrSuEdge* findFEdgeNot( const pmg::Vertex* vert ) const;
    void addFEdge(    const FrSuEdge* fEdge );
    void removeFEdge( const FrSuEdge* fEdge );
    bool isFEdgesFull() const;

    bool contains( const FrSuEdge* fEdge ) const;

    Facet(const Crystallite* relatedCrys, const pmg::Facet* facet);


private:
    Crystallite* m_relatedCrys;
};

} // namespace surface
} // namespace front
} // namespace pmg
