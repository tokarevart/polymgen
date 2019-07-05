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
    pmg::Face* x;
    // TODO: specify that filling must occur in method void setFEdges( ... );
    // TODO: add method to change 1 fedge like bool changeFEdge( const front::Edge* oldFEdge, const front::Edge* newFEdge ); in case of front splitting
    std::array<front::Edge*, 3> front_edges = { nullptr, nullptr, nullptr };
    vec3 normal;

    // TODO: add method std::array<pmg::Vert*, 3> verts() const;

    vec3 compute_normal();
    vec3 center();
    real_t quality();

    front::Vert* find_front_vert_not(const front::Edge* fedge) const;
    front::Edge* find_front_edge(const front::Vert* v0, const front::Vert* v1) const;
    front::Edge* find_front_edge(const pmg::Edge* edge) const;
    front::Edge* find_front_edge(const pmg::Vert* v0, const pmg::Vert* v1) const;
    front::Edge* find_front_edge_not(const front::Vert* fvert) const;
    front::Edge* find_front_edge_not(const pmg::Vert* vert) const;
    void add_front_edge(const front::Edge* fedge);
    void remove_front_edge(const front::Edge* fedge);
    bool front_edges_full() const;

    bool contains(const front::Edge* fedge) const;
    bool contains(const front::Vert* fvert) const;
    bool contains(const   pmg::Edge* edge) const;
    bool contains(const   pmg::Vert* vert) const;

    Face(const Polyhedron* related_polyhedron, const front::Edge* fedge0, const front::Edge* fedge1, const front::Edge* fedge2);
    Face(const Polyhedron* related_polyhedron, const pmg::Face* face);


private:
    // TODO: it's not good to store it
    Polyhedron* m_related_polyhedron;
};

} // namespace pmg::front
