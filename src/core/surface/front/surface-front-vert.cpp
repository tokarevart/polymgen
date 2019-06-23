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
using spt::vec3;


void front::Vert::refreshAngleData() {
    m_needAngleProcessing = true;
    m_needComplexityProcessing = true;
}


real_t front::Vert::complexity() {
    if (!m_needComplexityProcessing)
        return m_complexity;

    return computeComplexity();
}


real_t front::Vert::angle() {
    if (!m_needAngleProcessing)
        return m_angle;

    return computeAngle();
}


real_t front::Vert::computeComplexity() {
    m_needComplexityProcessing = false;
    auto adj_edges = findAdjEdges();
    real_t adj_edges_av_len = static_cast<real_t>(0.5) * (std::get<0>(adj_edges)->edge->magnitude()
                                                          + std::get<1>(adj_edges)->edge->magnitude());
    return m_complexity = m_relatedSurfaceFace->preferredLength() / adj_edges_av_len + K_ALPHA * PI / computeAngle();
}


real_t front::Vert::computeAngle() {
    auto adj_edges = findAdjEdges();
    real_t normals_cos = vec3::dot(std::get<0>(adj_edges)->normal, std::get<1>(adj_edges)->normal);

    m_needAngleProcessing = false;
    return m_angle = spt::algs::cpaTime(
        std::get<0>(adj_edges)->computeCenter(), std::get<0>(adj_edges)->normal,
        std::get<1>(adj_edges)->computeCenter(), std::get<1>(adj_edges)->normal) < static_cast<real_t>(1e-6) ?
        std::acos(normals_cos) + PI :
        std::acos(-normals_cos);
}




pair_ee front::Vert::findAdjEdges() const {
    pair_ee res(nullptr, nullptr);
    std::size_t i = 0;
    for (auto& fedge : m_relatedSurfaceFace->frontEdges()) {
        if (fedge->edge->contains(vert)) {
            if (i == 0) {
                res.first = fedge; i++;
            } else {
                res.second = fedge; return res;
            }
        }
    }

    throw std::logic_error("Function pmg::surface::front::Vert::findAdjEdges didn't find 2 adjacent front edges.");
}


pair_vv front::Vert::oppVerts() const {
    auto adj_edges = findAdjEdges();
    return pair_vv(adj_edges.first->edge->findNot(vert),
                   adj_edges.second->edge->findNot(vert));
}




front::Vert::Vert(const pmg::surface::Face* relatedSurfaceFace, const pmg::Vert* vert)
    : vert(const_cast<pmg::Vert*>(vert)), m_relatedSurfaceFace(const_cast<surface::Face*>(relatedSurfaceFace)) {}
