#include "spatial-objs/front/surface/front-surface-edge.h"
#include <cmath>
#include "helpers/spatial-algs/spatial-algs.h"


#define K_ALPHA 2.5


using FrSuFacet = pmg::front::surface::Facet;
using FrSuEdge  = pmg::front::surface::Edge;
using pair_vv = std::pair<pmg::Vertex*, pmg::Vertex*>;
using pair_ff = std::pair<FrSuFacet*, FrSuFacet*>;
using Vec   = tva::Vec;
using Point = tva::Point;




void FrSuEdge::refreshAngleData()
{
    m_needExCosProcessing      = true;
    m_needComplexityProcessing = true;
}


double FrSuEdge::complexity()
{
    if (!m_needComplexityProcessing)
        return m_complexity;

    return computeComplexity();
}


double FrSuEdge::angleExCos()
{
    if (!m_needExCosProcessing)
        return m_exCos;

    return computeAngleExCos();
}


double FrSuEdge::computeComplexity()
{
    m_needComplexityProcessing = false;
    return m_complexity = m_relatedCrys->preferredLength() / edge->magnitude() + K_ALPHA * M_PI / computeAngle();
}


double FrSuEdge::computeAngleExCos()
{
    auto adj_facets = getAdjFFacets();

    double normals_cos = tva::Vec::dot(std::get<0>(adj_facets)->normal, std::get<1>(adj_facets)->normal);

    m_needExCosProcessing = false;
    return m_exCos = tva::spatalgs::cpaTime(
            std::get<0>(adj_facets)->computeCenter(), std::get<0>(adj_facets)->normal,
            std::get<1>(adj_facets)->computeCenter(), std::get<1>(adj_facets)->normal) < 1e-6 ?
        -2.0 + normals_cos :
        -normals_cos;
}


double FrSuEdge::computeAngle()
{
    return angleExCos() < -1.0 ?
        acos(m_exCos + 2.0) + M_PI :
        acos(m_exCos);
}




pair_vv FrSuEdge::findOppVerts()
{
    auto adj_facets = getAdjFFacets();

    return { std::get<0>(adj_facets)->facet->findVertNot(edge),
             std::get<1>(adj_facets)->facet->findVertNot(edge) };
}


FrSuEdge* FrSuEdge::findOppEdge()
{
    auto opp_verts = findOppVerts();
    std::vector<FrSuEdge*> opp_fedges;
    for (auto& f_edge : m_relatedCrys->frontEdges())
    {
        if (((f_edge->edge->verts[0] == std::get<0>(opp_verts) &&
              f_edge->edge->verts[1] == std::get<1>(opp_verts)) ||
             (f_edge->edge->verts[1] == std::get<0>(opp_verts) &&
              f_edge->edge->verts[0] == std::get<1>(opp_verts))))
        {
            opp_fedges.push_back(f_edge);
        }
    }

    if (opp_fedges.size() == 1)
        return opp_fedges.front();

    std::vector<pair_ff> adj_ffacets_vec;
    adj_ffacets_vec.reserve(opp_fedges.size());
    for (auto& fedge : opp_fedges)
        adj_ffacets_vec.push_back(fedge->getAdjFFacets());

    Point main_vert_pos = edge->verts[0]->pos();
    Point main_vert_proj = tva::spatalgs::project(main_vert_pos, opp_verts.first->pos(), opp_verts.second->pos());
    Vec main_vec = main_vert_pos - main_vert_proj;

    FrSuEdge* max_cos_fedge = nullptr;
    double max_cos = -1.0;
    for (size_t i = 0; i < adj_ffacets_vec.size(); i++)
    {
        Point adj_opp_pos0 = adj_ffacets_vec[i].first->facet->findVertNot(opp_fedges.front()->edge)->pos();
        Point adj_opp_pos1 = adj_ffacets_vec[i].second->facet->findVertNot(opp_fedges.front()->edge)->pos();
        Point adj_opp_proj0 = tva::spatalgs::project(adj_opp_pos0, opp_verts.first->pos(), opp_verts.second->pos());
        Point adj_opp_proj1 = tva::spatalgs::project(adj_opp_pos1, opp_verts.first->pos(), opp_verts.second->pos());
        Vec adj_vec0 = adj_opp_pos0 - adj_opp_proj0;
        Vec adj_vec1 = adj_opp_pos1 - adj_opp_proj1;
        double cos0 = Vec::cos(adj_vec0, main_vec);
        double cos1 = Vec::cos(adj_vec1, main_vec);
        double max_cur_cos = std::max(cos0, cos1);
        if (max_cur_cos > max_cos)
        {
            max_cos = max_cur_cos;
            max_cos_fedge = opp_fedges[i];
        }
    }

    return max_cos_fedge;
}




pair_ff FrSuEdge::getAdjFFacets()
{
    if (!isAdjFacetsFull())
        fillAdjFFacets();

    return m_adjFFacets;
}




bool FrSuEdge::addAdjFFacet(const FrSuFacet* fFacet)
{
    if (!std::get<0>(m_adjFFacets))
        std::get<0>(m_adjFFacets)  = const_cast<FrSuFacet*>(fFacet);
    else if (!std::get<1>(m_adjFFacets))
        std::get<1>(m_adjFFacets) = const_cast<FrSuFacet*>(fFacet);
    else
        return false;

    return true;
}


bool FrSuEdge::removeAdjFFacet(const FrSuFacet* fFacet)
{
    if (std::get<0>(m_adjFFacets) == const_cast<FrSuFacet*>(fFacet))
        std::get<0>(m_adjFFacets) = nullptr;
    else if (std::get<1>(m_adjFFacets) == const_cast<FrSuFacet*>(fFacet))
        std::get<1>(m_adjFFacets) = nullptr;
    else
        return false;

    return true;
}


bool FrSuEdge::adjFFacetsContains(const FrSuFacet* fFacet) const
{
    return std::get<0>(m_adjFFacets) == const_cast<FrSuFacet*>(fFacet) ||
            std::get<1>(m_adjFFacets) == const_cast<FrSuFacet*>(fFacet);
}


void FrSuEdge::fillAdjFFacets(const FrSuFacet* fFacet0, const FrSuFacet* fFacet1)
{
    m_adjFFacets.first  = const_cast<FrSuFacet*>(fFacet0);
    m_adjFFacets.second = const_cast<FrSuFacet*>(fFacet1);
}




FrSuEdge::Edge(const Crystallite* relatedCrys, const pmg::Edge* edge)
    : edge(const_cast<pmg::Edge*>(edge)), m_relatedCrys(const_cast<Crystallite*>(relatedCrys)) {}




bool FrSuEdge::isAdjFacetsFull()
{
    return std::get<0>(m_adjFFacets) && std::get<1>(m_adjFFacets);
}


pair_ff FrSuEdge::fillAdjFFacets()
{
    for (auto& f_facet : m_relatedCrys->frontFacets())
    {
        if (f_facet->contains(this))
        {
            if (adjFFacetsContains(f_facet))
                continue;

            addAdjFFacet(f_facet);
            if (isAdjFacetsFull())
                break;
        }
    }

    if (!isAdjFacetsFull())
        throw std::logic_error("pmg::front::surface::Edge::fillAdjFFacets didn't find 2 adjacent front facets.");

    return m_adjFFacets;
}
