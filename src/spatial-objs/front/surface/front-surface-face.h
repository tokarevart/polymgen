// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/polyhedron.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace front {
namespace surface {

class Face
{
    using FrSuEdge = front::surface::Edge;

public:
    pmg::Face* face;
    FrSuEdge* fEdges[3] { nullptr, nullptr, nullptr };
    Vec normal;

    Vec    computeNormal();
    Vec    computeCenter();
    real_t computeQuality();

    FrSuEdge* findFEdge( const pmg::Edge* edge ) const;
    FrSuEdge* findFEdge( const pmg::Vertex* v0, const pmg::Vertex* v1 ) const;
    FrSuEdge* findFEdgeNot( const pmg::Vertex* vert ) const;
    void addFEdge(    const FrSuEdge* fEdge );
    void removeFEdge( const FrSuEdge* fEdge );
    bool isFEdgesFull() const;

    bool contains( const FrSuEdge* fEdge ) const;

    Face(const Polyhedron* relatedPolyhedron, const pmg::Face* face);


private:
    Polyhedron* m_relatedPolyhedron;
};

} // namespace surface
} // namespace front
} // namespace pmg
