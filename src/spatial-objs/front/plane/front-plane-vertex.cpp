// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/front/plane/front-plane-vertex.h"
#include <stdexcept>
#include "helpers/spatial-algs/spatial-algs.h"


#define K_ALPHA 4.0


using FrPlEdge   = pmg::front::plane::Edge;
using FrPlVertex = pmg::front::plane::Vertex;
using pair_ee  = std::pair<FrPlEdge*, FrPlEdge*>;
using pair_vv  = std::pair<pmg::Vertex*, pmg::Vertex*>;




void FrPlVertex::refreshAngleData()
{
    m_needExCosProcessing      = true;
    m_needComplexityProcessing = true;
}


double FrPlVertex::complexity()
{
    if (!m_needComplexityProcessing)
        return m_complexity;

    return computeComplexity();
}


double FrPlVertex::angleExCos()
{
    if (!m_needExCosProcessing)
        return m_exCos;

    return computeAngleExCos();
}


double FrPlVertex::computeComplexity()
{
    m_needComplexityProcessing = false;
    auto adj_edges = findAdjEdges();
    double adj_edges_av_len = 0.5 * (  std::get<0>(adj_edges)->edge->magnitude()
                                     + std::get<1>(adj_edges)->edge->magnitude());
    return m_complexity = m_relatedShellFacet->preferredLength() / adj_edges_av_len + K_ALPHA * M_PI / computeAngle();
}


double FrPlVertex::computeAngleExCos()
{
    auto adj_edges = findAdjEdges();

    double normals_cos = tva::Vec::dot(std::get<0>(adj_edges)->normal, std::get<1>(adj_edges)->normal);

    m_needExCosProcessing = false;
    return m_exCos = tva::spatalgs::cpaTime(
            std::get<0>(adj_edges)->computeCenter(), std::get<0>(adj_edges)->normal,
            std::get<1>(adj_edges)->computeCenter(), std::get<1>(adj_edges)->normal) < 1e-6 ?
        -2.0 + normals_cos :
        -normals_cos;
}


double FrPlVertex::computeAngle()
{
    return angleExCos() < -1.0 ?
        acos(m_exCos + 2.0) + M_PI :
                acos(m_exCos);
}




pair_ee FrPlVertex::findAdjEdges() const
{
    pair_ee res(nullptr, nullptr);
    int i = 0;
    for (auto& fedge : m_relatedShellFacet->frontEdges())
    {
        if (fedge->edge->contains(vert))
        {
            if (i == 0) { res.first  = fedge; i++; }
            else        { res.second = fedge; return res; }
        }
    }

    throw std::logic_error("Function pmg::front::plane::Vertex::findAdjEdges didn't find 2 adjacent front edges.");
}


pair_vv FrPlVertex::findOppVerts() const
{
    auto adj_edges = findAdjEdges();
    return pair_vv(adj_edges.first->edge->findNot(vert),
                   adj_edges.second->edge->findNot(vert));
}




FrPlVertex::Vertex(const shell::Facet* relatedShellFacet, const pmg::Vertex* vert)
    : vert(const_cast<pmg::Vertex*>(vert)), m_relatedShellFacet(const_cast<shell::Facet*>(relatedShellFacet)) {}
