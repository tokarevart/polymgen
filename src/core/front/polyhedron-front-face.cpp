// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "face.h"
#include <stdexcept>
#include "../../helpers/spatial/algs.h"

using namespace pmg;
using vec3 = spt::vec<3, real_t>;


vec3 front::Face::compute_normal() {
    vec3 center = this->center();
    vec3 third_pos = x->find_vert_not(x->edges[0])->pos();
    vec3 loc_normal = spt::cross(
        x->edges[0]->verts[0]->pos() - third_pos,
        x->edges[0]->verts[1]->pos() - third_pos).normalize();

    // I don't know why it led to better(correct) result.
    vec3 test_normal_correct_intersect = loc_normal + vec3(static_cast<real_t>(2.1632737147),
                                                           static_cast<real_t>(1.488313178),
                                                           static_cast<real_t>(-0.71123534278))
        * static_cast<real_t>(1e-3);

    std::size_t intersects_num = 0;
    for (auto& fface : m_related_polyhedron->front_faces()) {
        if (fface == this)
            continue;

        if (spt::does_ray_intersect_triangle(
            center, test_normal_correct_intersect,
            fface->x->edges[0]->verts[0]->pos(),
            fface->x->edges[0]->verts[1]->pos(),
            fface->x->find_vert_not(fface->x->edges[0])->pos()))
            intersects_num++;
    }

    if (loc_normal.sqr_magnitude() < static_cast<real_t>(1e-6))
        throw std::exception();

    return normal = intersects_num % 2 == 1 ? loc_normal : -loc_normal;
}


vec3 front::Face::center() {
    return x->center();
}


real_t front::Face::quality() {
    return x->quality();
}




front::Edge* front::Face::find_front_edge(const pmg::Edge* edge) const {
    for (auto& l_fedge : front_edges)
        if (l_fedge->x == edge)
            return l_fedge;

    return nullptr;
}


front::Edge* front::Face::find_front_edge(const pmg::Vert* v0, const pmg::Vert* v1) const {
    for (auto& l_fedge : front_edges)
        if ((l_fedge->x->verts[0] == v0 && l_fedge->x->verts[1] == v1) ||
            (l_fedge->x->verts[0] == v1 && l_fedge->x->verts[1] == v0))
            return l_fedge;

    return nullptr;
}


front::Edge* front::Face::find_front_edge_not(const pmg::Vert* vert) const {
    for (auto& l_fedge : front_edges)
        if (!l_fedge->x->contains(vert))
            return l_fedge;

    return nullptr;
}


void front::Face::add_front_edge(const front::Edge* fedge) {
    if (!front_edges[0])
        front_edges[0] = const_cast<front::Edge*>(fedge);
    else if (!front_edges[1])
        front_edges[1] = const_cast<front::Edge*>(fedge);
    else if (!front_edges[2])
        front_edges[2] = const_cast<front::Edge*>(fedge);
    else
        throw std::logic_error("pmg::front::Face::add_front_edge can't add fedge when front_edges already full.");
}


void front::Face::remove_front_edge(const front::Edge* fedge) {
    for (auto& l_fedge : front_edges) {
        if (l_fedge == fedge) {
            l_fedge = nullptr;
            return;
        }
    }

    throw std::logic_error("pmg::front::Face::remove_front_edge can't remove fedge because there is no one.");
}


bool front::Face::front_edges_full() const {
    if (front_edges[0] && front_edges[1] && front_edges[2])
        return true;

    return false;
}


bool front::Face::contains(const front::Edge* fedge) const {
    if (front_edges[0] == fedge ||
        front_edges[1] == fedge ||
        front_edges[2] == fedge)
        return true;

    return false;
}


bool front::Face::contains(const pmg::Edge* edge) const {
    return x->contains(edge);
}


bool front::Face::contains(const pmg::Vert* vert) const {
    return x->contains(vert);
}




front::Face::Face(const Polyhedron* related_polyhedron, const pmg::Face* face)
    : x(const_cast<pmg::Face*>(face)), m_related_polyhedron(const_cast<Polyhedron*>(related_polyhedron)) {}
