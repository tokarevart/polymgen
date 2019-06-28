// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "../polyhedron.h"
#include "edge.h"
#include "../../helpers/spatial/vec.h"
#include "../../real-type.h"

#include "../../definitions.h"


namespace pmg::front {

class Face {
    using vec3 = spt::vec<3, real_t>;

public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    pmg::Face* face;
    // TODO: specify that filling must occur in method void setFEdges( ... );
    // TODO: add method to change 1 front_edge like bool changeFEdge( const front::Edge* oldFEdge, const front::Edge* newFEdge ); in case of front splitting
    std::array<front::Edge*, 3> front_edges = { nullptr, nullptr, nullptr };
    vec3 normal;

    // TODO: add method std::array<pmg::Vert*, 3> verts() const;

    vec3 compute_normal();
    vec3 center();
    real_t quality();

    front::Edge* findFEdge(const pmg::Edge* edge) const;
    front::Edge* findFEdge(const pmg::Vert* v0, const pmg::Vert* v1) const;
    front::Edge* findFEdgeNot(const pmg::Vert* vert) const;
    void addFEdge(const front::Edge* front_edge);
    void removeFEdge(const front::Edge* front_edge);
    bool fEdgesFull() const;

    bool contains(const front::Edge* front_edge) const;
    bool contains(const   pmg::Edge* edge) const;
    bool contains(const   pmg::Vert* vert) const;

    Face(const Polyhedron* related_polyhedron, const pmg::Face* face);


private:
    // TODO: it's not good to store it
    Polyhedron* m_related_polyhedron;
};

} // namespace pmg::front
