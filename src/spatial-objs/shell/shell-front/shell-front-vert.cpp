// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/shell/shell-front/shell-front-vert.h"
#include <stdexcept>
#include "helpers/spatial-algs/spatial-algs.h"

#define PI static_cast<real_t>(M_PI)

#define K_ALPHA static_cast<real_t>(4.0)


using FrPlEdge = pmg::shell::front::Edge;
using FrPlVert = pmg::shell::front::Vert;
using pair_ee  = std::pair<FrPlEdge*, FrPlEdge*>;
using pair_vv  = std::pair<pmg::Vert*, pmg::Vert*>;




void FrPlVert::refreshAngleData()
{
    m_needAngleProcessing      = true;
    m_needComplexityProcessing = true;
}


real_t FrPlVert::complexity()
{
    if (!m_needComplexityProcessing)
        return m_complexity;

    return computeComplexity();
}


real_t FrPlVert::angle()
{
    if (!m_needAngleProcessing)
        return m_angle;

    return computeAngle();
}


real_t FrPlVert::computeComplexity()
{
    m_needComplexityProcessing = false;
    auto adj_edges = findAdjEdges();
    real_t adj_edges_av_len = static_cast<real_t>(0.5) * (  std::get<0>(adj_edges)->edge->magnitude()
                                     + std::get<1>(adj_edges)->edge->magnitude());
    return m_complexity = m_relatedShellFace->preferredLength() / adj_edges_av_len + K_ALPHA * PI / computeAngle();
}


real_t FrPlVert::computeAngle()
{
    auto adj_edges = findAdjEdges();
    real_t normals_cos = vec3::dot(std::get<0>(adj_edges)->normal, std::get<1>(adj_edges)->normal);

    m_needAngleProcessing = false;
    return m_angle = spatalgs::cpaTime(
            std::get<0>(adj_edges)->computeCenter(), std::get<0>(adj_edges)->normal,
            std::get<1>(adj_edges)->computeCenter(), std::get<1>(adj_edges)->normal) < static_cast<real_t>(1e-6) ?
        std::acos(normals_cos) + PI :
        std::acos(-normals_cos);
}




pair_ee FrPlVert::findAdjEdges() const
{
    pair_ee res(nullptr, nullptr);
    int i = 0;
    for (auto& fedge : m_relatedShellFace->frontEdges())
    {
        if (fedge->edge->contains(vert))
        {
            if (i == 0) { res.first  = fedge; i++; }
            else        { res.second = fedge; return res; }
        }
    }

    throw std::logic_error("Function pmg::shell::front::Vert::findAdjEdges didn't find 2 adjacent front edges.");
}


pair_vv FrPlVert::oppVerts() const
{
    auto adj_edges = findAdjEdges();
    return pair_vv(adj_edges.first->edge->findNot(vert),
                   adj_edges.second->edge->findNot(vert));
}




FrPlVert::Vert(const shell::Face* relatedShellFace, const pmg::Vert* vert)
    : vert(const_cast<pmg::Vert*>(vert)), m_relatedShellFace(const_cast<shell::Face*>(relatedShellFace)) {}
