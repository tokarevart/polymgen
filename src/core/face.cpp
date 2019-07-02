// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "face.h"
#include "../helpers/spatial/algs.h"

#define ONE_3 static_cast<real_t>(0.3333333333333333)

using vec3 = spt::vec<3, real_t>;


vec3 pmg::Face::center() const {
    return (edges[0]->verts[0]->pos() +
            edges[0]->verts[1]->pos() +
            find_vert_not(edges[0])->pos()) * ONE_3;
}


real_t pmg::Face::quality() const {
    return shortest_edge()->sqr_magnitude() / longest_edge()->sqr_magnitude();
}


real_t pmg::Face::area() const {
    vec3 vec0 = edges[0]->verts[1]->pos() - edges[0]->verts[0]->pos();
    vec3 vec1 = edges[1]->verts[1]->pos() - edges[1]->verts[0]->pos();
    return spt::cross(vec0, vec1).magnitude() * static_cast<real_t>(0.5);
}


pmg::Vert* pmg::Face::find_vert_not(const pmg::Edge* edge) const {
    for (auto& face_edge : edges) {
        if (edge != face_edge) {
            if (!edge->contains(face_edge->verts[0]))
                return face_edge->verts[0];
            else
                return face_edge->verts[1];
        }
    }

    return nullptr;
}


pmg::Edge* pmg::Face::find_edge_not(const pmg::Vert* vert) const {
    for (auto& edge : edges)
        if (!edge->contains(vert))
            return edge;

    return nullptr;
}


pmg::Edge* pmg::Face::find_edge(const pmg::Vert* vert0, const pmg::Vert* vert1) const {
    for (auto& edge : edges)
        if ((edge->verts[0] == vert0 &&
             edge->verts[1] == vert1) ||
             (edge->verts[0] == vert1 &&
              edge->verts[1] == vert0))
            return edge;

    return nullptr;
}


pmg::Edge* pmg::Face::shortest_edge() const {
    if (edges[0]->sqr_magnitude() < edges[1]->sqr_magnitude()) {
        if (edges[0]->sqr_magnitude() < edges[2]->sqr_magnitude())
            return edges[0];
        else
            return edges[2];
    } else {
        if (edges[1]->sqr_magnitude() < edges[2]->sqr_magnitude())
            return edges[1];
        else
            return edges[2];
    }
}


pmg::Edge* pmg::Face::longest_edge() const {
    if (edges[0]->sqr_magnitude() > edges[1]->sqr_magnitude()) {
        if (edges[0]->sqr_magnitude() > edges[2]->sqr_magnitude())
            return edges[0];
        else
            return edges[2];
    } else {
        if (edges[1]->sqr_magnitude() > edges[2]->sqr_magnitude())
            return edges[1];
        else
            return edges[2];
    }
}


bool pmg::Face::contains(const pmg::Edge* edge) const {
    for (auto& edge0 : edges)
        if (edge0 == edge)
            return true;

    return false;
}


bool pmg::Face::contains(const pmg::Vert* vert) const {
    for (auto& edge : edges)
        if (edge->contains(vert))
            return true;

    return false;
}


pmg::Face::Face(const pmg::Edge* edge0, const pmg::Edge* edge1, const pmg::Edge* edge2) {
    edges[0] = const_cast<pmg::Edge*>(edge0);
    edges[1] = const_cast<pmg::Edge*>(edge1);
    edges[2] = const_cast<pmg::Edge*>(edge2);
}


pmg::Face::Face(const pmg::Vert* vert0, const pmg::Vert* vert1, const pmg::Vert* vert2) {
    edges[0] = new pmg::Edge(vert0, vert1);
    edges[1] = new pmg::Edge(vert1, vert2);
    edges[2] = new pmg::Edge(vert2, vert0);
}
