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
    return m_complexity = m_related_polyhedron->preferred_length() / edge->magnitude() + K_ALPHA * PI / angle();
}


real_t front::Edge::compute_angle() {
    auto adj_faces = adj_ffaces();
    real_t normals_cos = spt::dot(std::get<0>(adj_faces)->normal, std::get<1>(adj_faces)->normal);

    m_need_angle_processing = false;
    return m_angle = spt::cpa_time(
        std::get<0>(adj_faces)->center(), std::get<0>(adj_faces)->normal,
        std::get<1>(adj_faces)->center(), std::get<1>(adj_faces)->normal) < static_cast<real_t>(1e-6) ?
        std::acos(normals_cos) + PI :
        std::acos(-normals_cos);
}


pair_vv front::Edge::opp_verts() {
    auto adj_faces = adj_ffaces();

    return { std::get<0>(adj_faces)->face->find_vert_not(edge),
             std::get<1>(adj_faces)->face->find_vert_not(edge) };
}


FrSuEdge* front::Edge::find_opp_edge() {
    auto opp_verts_l = opp_verts();
    std::vector<front::Edge*> opp_fedges;
    for (auto& f_edge : m_related_polyhedron->front_edges()) {
        if (((f_edge->edge->verts[0] == std::get<0>(opp_verts_l) &&
              f_edge->edge->verts[1] == std::get<1>(opp_verts_l)) ||
              (f_edge->edge->verts[1] == std::get<0>(opp_verts_l) &&
               f_edge->edge->verts[0] == std::get<1>(opp_verts_l)))) {
            opp_fedges.push_back(f_edge);
        }
    }

    if (opp_fedges.size() == 1)
        return opp_fedges.front();

    std::vector<pair_ff> adj_ffaces_vec;
    adj_ffaces_vec.reserve(opp_fedges.size());
    for (auto& f_edge : opp_fedges)
        adj_ffaces_vec.push_back(f_edge->adj_ffaces());

    vec3 main_vert_pos = edge->verts[0]->pos();
    vec3 main_vert_proj = spt::project(main_vert_pos, opp_verts_l.first->pos(), opp_verts_l.second->pos());
    vec3 main_vec = main_vert_pos - main_vert_proj;

    front::Edge* max_cos_fedge = nullptr;
    real_t max_cos = -1.0;
    for (std::size_t i = 0; i < adj_ffaces_vec.size(); i++) {
        vec3 adj_opp_pos0 = adj_ffaces_vec[i].first->face->find_vert_not(opp_fedges.front()->edge)->pos();
        vec3 adj_opp_pos1 = adj_ffaces_vec[i].second->face->find_vert_not(opp_fedges.front()->edge)->pos();
        vec3 adj_opp_proj0 = spt::project(adj_opp_pos0, opp_verts_l.first->pos(), opp_verts_l.second->pos());
        vec3 adj_opp_proj1 = spt::project(adj_opp_pos1, opp_verts_l.first->pos(), opp_verts_l.second->pos());
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


bool front::Edge::add_adj_fface(const front::Face* fFace) {
    if (!std::get<0>(m_adj_ffaces))
        std::get<0>(m_adj_ffaces) = const_cast<front::Face*>(fFace);
    else if (!std::get<1>(m_adj_ffaces))
        std::get<1>(m_adj_ffaces) = const_cast<front::Face*>(fFace);
    else
        return false;

    return true;
}


bool front::Edge::remove_adj_fface(const front::Face* fFace) {
    if (std::get<0>(m_adj_ffaces) == const_cast<front::Face*>(fFace))
        std::get<0>(m_adj_ffaces) = nullptr;
    else if (std::get<1>(m_adj_ffaces) == const_cast<front::Face*>(fFace))
        std::get<1>(m_adj_ffaces) = nullptr;
    else
        return false;

    return true;
}


bool front::Edge::adj_ffaces_contains(const front::Face* fFace) const {
    return std::get<0>(m_adj_ffaces) == const_cast<front::Face*>(fFace) ||
        std::get<1>(m_adj_ffaces) == const_cast<front::Face*>(fFace);
}


void front::Edge::fill_adj_ffaces(const front::Face* fFace0, const front::Face* fFace1) {
    m_adj_ffaces.first = const_cast<front::Face*>(fFace0);
    m_adj_ffaces.second = const_cast<front::Face*>(fFace1);
}




front::Edge::Edge(const Polyhedron* related_polyhedron, const pmg::Edge* edge)
    : edge(const_cast<pmg::Edge*>(edge)), m_related_polyhedron(const_cast<Polyhedron*>(related_polyhedron)) {}




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
