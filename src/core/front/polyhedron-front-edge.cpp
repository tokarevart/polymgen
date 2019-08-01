// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "edge.h"
#include <cmath>
#include "../../helpers/mathconsts.h"
#include "../../helpers/spatial/algs.h"


#define K_ALPHA static_cast<real_t>(4.0)


using namespace pmg;
using FrSuEdge = pmg::front::Edge;
using pair_vv = std::pair<pmg::Vert*, pmg::Vert*>;
using pair_ff = std::pair<front::Face*, front::Face*>;
using vec3 = spt::vec<3, real_t>;


void pmg::front::Edge::add_front_vert(const front::Vert* fvert) {
    if(!front_verts[0])
        front_verts[0] = const_cast<front::Vert*>(fvert);
    else if (!front_verts[1])
        front_verts[1] = const_cast<front::Vert*>(fvert);
    else
        throw std::logic_error("pmg::front::Edge::add_front_vert can't add fvert when front_edges already full.");
}


void pmg::front::Edge::remove_front_vert(const front::Vert* fvert) {
    for (auto& l_fvert : front_verts)
        if (l_fvert == fvert) {
            l_fvert = nullptr;
            return;
        }

    throw std::logic_error("pmg::front::Edge::remove_front_vert can't remove fvert because there is no one.");
}


bool pmg::front::Edge::front_verts_full() const {
    if (front_verts[0] && front_verts[1])
        return true;

    return false;
}


void front::Edge::refresh_angle_data() {
    m_need_angle_processing = true;
    m_need_complexity_processing = true;
}


real_t front::Edge::complexity() {
    if (!m_need_complexity_processing)
        return m_complexity;

    return compute_complexity();
}


real_t front::Edge::angle() {
    if (!m_need_angle_processing)
        return m_angle;

    return compute_angle();
}


real_t front::Edge::compute_complexity() {
    m_need_complexity_processing = false;
    return m_complexity = m_related_polyhedron->preferred_length() / x->magnitude() + K_ALPHA * mathconsts::PI / angle();
}


real_t front::Edge::compute_angle() {
    auto adj_faces = adj_ffaces();
    real_t normals_cos = spt::dot(std::get<0>(adj_faces)->normal, std::get<1>(adj_faces)->normal);

    m_need_angle_processing = false;
    return m_angle = spt::cpa_time(
        std::get<0>(adj_faces)->center(), std::get<0>(adj_faces)->normal,
        std::get<1>(adj_faces)->center(), std::get<1>(adj_faces)->normal) < static_cast<real_t>(1e-6) ?
        std::acos(normals_cos) + mathconsts::PI :
        std::acos(-normals_cos);
}


pair_vv front::Edge::opp_verts() {
    auto adj_faces = adj_ffaces();

    return { std::get<0>(adj_faces)->x->find_vert_not(x),
             std::get<1>(adj_faces)->x->find_vert_not(x) };
}


FrSuEdge* front::Edge::find_opp_edge() {
    auto l_opp_verts = opp_verts();
    std::vector<front::Edge*> opp_fedges;
    for (auto& l_fedge : m_related_polyhedron->front_edges()) {
        if (((l_fedge->x->verts[0] == std::get<0>(l_opp_verts) &&
              l_fedge->x->verts[1] == std::get<1>(l_opp_verts)) ||
              (l_fedge->x->verts[1] == std::get<0>(l_opp_verts) &&
               l_fedge->x->verts[0] == std::get<1>(l_opp_verts)))) {
            opp_fedges.push_back(l_fedge);
        }
    }

    if (opp_fedges.size() == 1)
        return opp_fedges.front();

    std::vector<pair_ff> adj_ffaces_vec;
    adj_ffaces_vec.reserve(opp_fedges.size());
    for (auto& l_fedge : opp_fedges)
        adj_ffaces_vec.push_back(l_fedge->adj_ffaces());

    vec3 main_vert_pos = x->verts[0]->pos();
    vec3 main_vert_proj = spt::project(main_vert_pos, l_opp_verts.first->pos(), l_opp_verts.second->pos());
    vec3 main_vec = main_vert_pos - main_vert_proj;

    front::Edge* max_cos_fedge = nullptr;
    real_t max_cos = -1.0;
    for (std::size_t i = 0; i < adj_ffaces_vec.size(); i++) {
        vec3 adj_opp_pos0 = adj_ffaces_vec[i].first->x->find_vert_not(opp_fedges.front()->x)->pos();
        vec3 adj_opp_pos1 = adj_ffaces_vec[i].second->x->find_vert_not(opp_fedges.front()->x)->pos();
        vec3 adj_opp_proj0 = spt::project(adj_opp_pos0, l_opp_verts.first->pos(), l_opp_verts.second->pos());
        vec3 adj_opp_proj1 = spt::project(adj_opp_pos1, l_opp_verts.first->pos(), l_opp_verts.second->pos());
        vec3 adj_vec0 = adj_opp_pos0 - adj_opp_proj0;
        vec3 adj_vec1 = adj_opp_pos1 - adj_opp_proj1;
        real_t cos0 = spt::cos(adj_vec0, main_vec);
        real_t cos1 = spt::cos(adj_vec1, main_vec);
        real_t max_cur_cos = std::max(cos0, cos1);
        if (max_cur_cos > max_cos) {
            max_cos = max_cur_cos;
            max_cos_fedge = opp_fedges[i];
        }
    }

    return max_cos_fedge;
}


pair_ff front::Edge::adj_ffaces() {
    if (!adj_faces_full())
        fill_adj_ffaces();

    return m_adj_ffaces;
}


bool front::Edge::add_adj_fface(const front::Face* fface) {
    if (!std::get<0>(m_adj_ffaces))
        std::get<0>(m_adj_ffaces) = const_cast<front::Face*>(fface);
    else if (!std::get<1>(m_adj_ffaces))
        std::get<1>(m_adj_ffaces) = const_cast<front::Face*>(fface);
    else
        return false;

    return true;
}


bool front::Edge::remove_adj_fface(const front::Face* fface) {
    if (std::get<0>(m_adj_ffaces) == const_cast<front::Face*>(fface))
        std::get<0>(m_adj_ffaces) = nullptr;
    else if (std::get<1>(m_adj_ffaces) == const_cast<front::Face*>(fface))
        std::get<1>(m_adj_ffaces) = nullptr;
    else
        return false;

    return true;
}


void pmg::front::Edge::clear_adj_ffaces() {
    std::get<0>(m_adj_ffaces) = nullptr;
    std::get<1>(m_adj_ffaces) = nullptr;
}


bool front::Edge::adj_ffaces_contains(const front::Face* fface) const {
    return std::get<0>(m_adj_ffaces) == const_cast<front::Face*>(fface) ||
        std::get<1>(m_adj_ffaces) == const_cast<front::Face*>(fface);
}


void front::Edge::fill_adj_ffaces(const front::Face* fFace0, const front::Face* fFace1) {
    m_adj_ffaces.first = const_cast<front::Face*>(fFace0);
    m_adj_ffaces.second = const_cast<front::Face*>(fFace1);
}




bool pmg::front::Edge::contains(const front::Vert* fvert) const {
    return front_verts[0] == fvert || front_verts[1] == fvert;
}

pmg::front::Edge::Edge(const Polyhedron* related_polyhedron, const front::Vert* fvert0, const front::Vert* fvert1)
    : m_related_polyhedron{ const_cast<Polyhedron*>(related_polyhedron) },
      front_verts{ const_cast<front::Vert*>(fvert0), const_cast<front::Vert*>(fvert1) },
      x{ new pmg::Edge(fvert0->x, fvert1->x) } {}


front::Edge::Edge(const Polyhedron* related_polyhedron, const pmg::Edge* edge)
    : x{ const_cast<pmg::Edge*>(edge) }, m_related_polyhedron{ const_cast<Polyhedron*>(related_polyhedron) } {}




bool front::Edge::adj_faces_full() {
    return std::get<0>(m_adj_ffaces) && std::get<1>(m_adj_ffaces);
}


pair_ff front::Edge::fill_adj_ffaces() {
    for (auto& fface : m_related_polyhedron->front_faces()) {
        if (fface->contains(this)) {
            if (adj_ffaces_contains(fface))
                continue;

            add_adj_fface(fface);
            if (adj_faces_full())
                break;
        }
    }

    if (!adj_faces_full())
        throw std::logic_error("pmg::front::Edge::fill_adj_ffaces didn't find 2 adjacent front Faces.");

    return m_adj_ffaces;
}
