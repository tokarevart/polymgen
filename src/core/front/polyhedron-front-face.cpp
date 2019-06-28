// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "face.h"
#include <stdexcept>
#include "../../helpers/spatial/algs.h"

using namespace pmg;
using vec3 = spt::vec<3, real_t>;


vec3 front::Face::compute_normal() {
    vec3 center = this->center();
    vec3 third_pos = face->find_vert_not(face->edges[0])->pos();
    vec3 loc_normal = spt::cross(
        face->edges[0]->verts[0]->pos() - third_pos,
        face->edges[0]->verts[1]->pos() - third_pos).normalize();

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
            fface->face->edges[0]->verts[0]->pos(),
            fface->face->edges[0]->verts[1]->pos(),
            fface->face->find_vert_not(fface->face->edges[0])->pos()))
            intersects_num++;
    }

    if (loc_normal.sqr_magnitude() < static_cast<real_t>(1e-6))
        throw std::exception();

    return normal = intersects_num % 2 == 1 ? loc_normal : -loc_normal;
}


vec3 front::Face::center() {
    return face->center();
}


real_t front::Face::quality() {
    return face->quality();
}




front::Edge* front::Face::findFEdge(const pmg::Edge* edge) const {
    for (auto& f_edge : front_edges)
        if (f_edge->edge == edge)
            return f_edge;

    return nullptr;
}


front::Edge* front::Face::findFEdge(const pmg::Vert* v0, const pmg::Vert* v1) const {
    for (auto& f_edge : front_edges)
        if ((f_edge->edge->verts[0] == v0 && f_edge->edge->verts[1] == v1) ||
            (f_edge->edge->verts[0] == v1 && f_edge->edge->verts[1] == v0))
            return f_edge;

    return nullptr;
}


front::Edge* front::Face::findFEdgeNot(const pmg::Vert* vert) const {
    for (auto& f_edge : front_edges)
        if (!f_edge->edge->contains(vert))
            return f_edge;

    return nullptr;
}


void front::Face::addFEdge(const front::Edge* front_edge) {
    if (!front_edges[0])
        front_edges[0] = const_cast<front::Edge*>(front_edge);
    else if (!front_edges[1])
        front_edges[1] = const_cast<front::Edge*>(front_edge);
    else if (!front_edges[2])
        front_edges[2] = const_cast<front::Edge*>(front_edge);
    else
        throw std::logic_error("pmg::front::Face::addFEdge can't add front_edge when front_edges already full.");
}


void front::Face::removeFEdge(const front::Edge* front_edge) {
    for (auto& f_edge : front_edges) {
        if (f_edge == front_edge) {
            f_edge = nullptr;
            return;
        }
    }

    throw std::logic_error("pmg::front::Face::removeFEdge can't remove front_edge because there is no one.");
}


bool front::Face::fEdgesFull() const {
    if (front_edges[0] && front_edges[1] && front_edges[2])
        return true;

    return false;
}


bool front::Face::contains(const front::Edge* front_edge) const {
    if (front_edges[0] == front_edge ||
        front_edges[1] == front_edge ||
        front_edges[2] == front_edge)
        return true;

    return false;
}


bool front::Face::contains(const pmg::Edge* edge) const {
    return face->contains(edge);
}


bool front::Face::contains(const pmg::Vert* vert) const {
    return face->contains(vert);
}




front::Face::Face(const Polyhedron* related_polyhedron, const pmg::Face* face)
    : face(const_cast<pmg::Face*>(face)), m_related_polyhedron(const_cast<Polyhedron*>(related_polyhedron)) {}
