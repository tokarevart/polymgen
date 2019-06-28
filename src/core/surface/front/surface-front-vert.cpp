// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "vert.h"
#include <stdexcept>
#include "../../../helpers/mathconsts.h"
#include "../../../helpers/spatial/algs.h"


#define K_ALPHA static_cast<real_t>(4.0)


using namespace pmg::surface;
using pair_ee  = std::pair<front::Edge*, front::Edge*>;
using pair_vv  = std::pair<pmg::Vert*, pmg::Vert*>;
using vec3 = spt::vec<3, real_t>;


void front::Vert::refresh_angle_data() {
    m_need_angle_processing = true;
    m_need_complexity_processing = true;
}


real_t front::Vert::complexity() {
    if (!m_need_complexity_processing)
        return m_complexity;

    return compute_complexity();
}


real_t front::Vert::angle() {
    if (!m_need_angle_processing)
        return m_angle;

    return compute_angle();
}


real_t front::Vert::compute_complexity() {
    m_need_complexity_processing = false;
    auto adj_edges = findAdjEdges();
    real_t adj_edges_av_len = static_cast<real_t>(0.5) * (std::get<0>(adj_edges)->edge->magnitude()
                                                          + std::get<1>(adj_edges)->edge->magnitude());
    return m_complexity = m_related_surface_face->preferred_length() / adj_edges_av_len + K_ALPHA * PI / compute_angle();
}


real_t front::Vert::compute_angle() {
    auto adj_edges = findAdjEdges();
    real_t normals_cos = spt::dot(std::get<0>(adj_edges)->normal, std::get<1>(adj_edges)->normal);

    m_need_angle_processing = false;
    return m_angle = spt::cpa_time(
        std::get<0>(adj_edges)->center(), std::get<0>(adj_edges)->normal,
        std::get<1>(adj_edges)->center(), std::get<1>(adj_edges)->normal) < static_cast<real_t>(1e-6) ?
        std::acos(normals_cos) + PI :
        std::acos(-normals_cos);
}




pair_ee front::Vert::findAdjEdges() const {
    pair_ee res(nullptr, nullptr);
    std::size_t i = 0;
    for (auto& f_edge : m_related_surface_face->front_edges()) {
        if (f_edge->edge->contains(vert)) {
            if (i == 0) {
                res.first = f_edge; i++;
            } else {
                res.second = f_edge; return res;
            }
        }
    }

    throw std::logic_error("Function pmg::surface::front::Vert::findAdjEdges didn't find 2 adjacent front edges.");
}


pair_vv front::Vert::opp_verts() const {
    auto adj_edges = findAdjEdges();
    return pair_vv(adj_edges.first->edge->findNot(vert),
                   adj_edges.second->edge->findNot(vert));
}




front::Vert::Vert(const pmg::surface::Face* relatedSurfaceFace, const pmg::Vert* vert)
    : vert(const_cast<pmg::Vert*>(vert)), m_related_surface_face(const_cast<surface::Face*>(relatedSurfaceFace)) {}
