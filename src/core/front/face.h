// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "core/polyhedron.h"
#include "core/front/edge.h"
#include "helpers/spatial/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg::front {

class Face
{
public:
    pmg::Face* face;
    // TODO: specify that filling must occur in method void assignFEdges( ... );
    // TODO: add method to change 1 fEdge like bool changeFEdge( const front::Edge* oldFEdge, const front::Edge* newFEdge ); in case of front splitting
    std::array<front::Edge*, 3> fEdges = { nullptr, };
    spt::vec3 normal;

    // TODO: add method std::array<pmg::Vert*, 3> verts() const;

    spt::vec3 computeNormal();
    spt::vec3 computeCenter();
    real_t computeQuality();

    front::Edge* findFEdge( const pmg::Edge* edge ) const;
    front::Edge* findFEdge( const pmg::Vert* v0, const pmg::Vert* v1 ) const;
    front::Edge* findFEdgeNot( const pmg::Vert* vert ) const;
    void addFEdge(    const front::Edge* fEdge );
    void removeFEdge( const front::Edge* fEdge );
    bool fEdgesFull() const;

    bool contains( const front::Edge* fEdge ) const;
    bool contains( const   pmg::Edge* edge  ) const;
    bool contains( const   pmg::Vert* vert  ) const;

    Face( const Polyhedron* relatedPolyhedron, const pmg::Face* face );


private:
    // TODO: it's not good to store it
    Polyhedron* m_relatedPolyhedron;
};

} // namespace pmg::front
