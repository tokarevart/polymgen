// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/polyhedron.h"
#include "spatial-objs/polyhedron-front/polyhedron-front-edge.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace front {

class Face
{
public:
    pmg::Face* face;
    std::array<front::Edge*, 3> fEdges = { nullptr, };
    vec3 normal;

    vec3   computeNormal();
    vec3   computeCenter();
    real_t computeQuality();

    front::Edge* findFEdge( const pmg::Edge* edge ) const;
    front::Edge* findFEdge( const pmg::Vert* v0, const pmg::Vert* v1 ) const;
    front::Edge* findFEdgeNot( const pmg::Vert* vert ) const;
    void addFEdge(    const front::Edge* fEdge );
    void removeFEdge( const front::Edge* fEdge );
    bool isFEdgesFull() const;

    bool contains( const   pmg::Vert* vert )  const;
    bool contains( const   pmg::Edge* edge )  const;
    bool contains( const front::Edge* fEdge ) const;

    Face(const Polyhedron* relatedPolyhedron, const pmg::Face* face);


private:
    Polyhedron* m_relatedPolyhedron;
};

} // namespace front
} // namespace pmg
