#include "Edge3.h"
#include <algorithm>
#include <iostream>
#include "helpers/spatialalgs/SpatialAlgs.h"


#define PI 3.141592653589793

#define EPS 1e-10
#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0[0], corner1[0], point[0]) && \
         BETWEEN(corner0[1], corner1[1], point[1]))


#define K_ALPHA 2.5




double Edge3::magnitude() const
{
    return (*verts[1] - *verts[0]).magnitude();
}

double Edge3::sqrMagnitude() const
{
    return (*verts[1] - *verts[0]).sqrMagnitude();
}




void Edge3::flip(std::list<Edge3*>& edgesList, std::list<Facet3*>& facetsList)
{
    Vertex3* around_nodes[2];
    Facet3* around_facets[2];
    find2AdjFacets(facetsList, around_facets[0], around_facets[1]);
    around_nodes[0] = around_facets[0]->findVertNotIncludedInEdge(this);
    around_nodes[1] = around_facets[1]->findVertNotIncludedInEdge(this);
    Edge3* new_edge = new Edge3(around_nodes[0], around_nodes[1]);
    edgesList.push_back(new_edge);

    facetsList.push_back(
        new Facet3(
            around_facets[0]->findEdge(around_nodes[0], verts[0]),
            around_facets[1]->findEdge(around_nodes[1], verts[0]),
            new_edge));
    facetsList.push_back(
        new Facet3(
            around_facets[0]->findEdge(around_nodes[0], verts[1]),
            around_facets[1]->findEdge(around_nodes[1], verts[1]),
            new_edge));

    facetsList.erase(std::find(facetsList.begin(), facetsList.end(), around_facets[0]));
    facetsList.erase(std::find(facetsList.begin(), facetsList.end(), around_facets[1]));
    delete around_facets[0];
    delete around_facets[1];

    edgesList.erase(std::find(edgesList.begin(), edgesList.end(), this));
    delete this;
}

bool Edge3::flipIfNeeded(std::list<Edge3*>& edgesList, std::list<Facet3*>& facetsList)
{
    if (needToFlip(facetsList))
    {
        flip(edgesList, facetsList);
        return true;
    }

    return false;
}

void Edge3::findAdjFacets(const std::list<Facet3*>& facetsList, std::list<Facet3*>& adjFacets)
{
    for (auto& facet : facetsList)
    {
        if (facet->contains(this))
            adjFacets.push_back(facet);
    }
}


void Edge3::find2AdjFacets(const std::list<Facet3*>& facetsList, Facet3*& facet0, Facet3*& facet1)
{
    bool not_found_yet = true;
    for (auto& facet : facetsList)
    {
        if (facet->contains(this))
        {
            if (not_found_yet)
            {
                facet0 = facet;
                not_found_yet = false;
            }
            else
            {
                facet1 = facet;
                return;
            }
        }
    }
}


void Edge3::make2Instead(std::list<Facet3*>& facetsList, std::list<Edge3*>& edgesList, std::list<Vertex3*>& vertsList)
{
    Vertex3* inner_vert = new Vertex3(verts[0]->getPos() + 0.5 * (*verts[1] - *verts[0]));
    vertsList.push_back(inner_vert);

    std::list<Facet3*> old_facets;
    findAdjFacets(facetsList, old_facets);

    Edge3* edge_halfs[2];
    edge_halfs[0] = new Edge3(verts[0], inner_vert);
    edge_halfs[1] = new Edge3(inner_vert, verts[1]);

    Vertex3* not_belongs_to_edge_vert;
    Edge3* dividing_edge;
    for (auto &facet : old_facets)
    {
        not_belongs_to_edge_vert = facet->findVertNotIncludedInEdge(this);
        dividing_edge = new Edge3(inner_vert, not_belongs_to_edge_vert);
        edgesList.push_back(dividing_edge);

        facetsList.push_back(
            new Facet3(
                facet->findEdge(verts[0], not_belongs_to_edge_vert),
                edge_halfs[0],
                dividing_edge));

        facetsList.push_back(
            new Facet3(
                facet->findEdge(verts[1], not_belongs_to_edge_vert),
                edge_halfs[1],
                dividing_edge));

        facetsList.erase(std::find(facetsList.begin(), facetsList.end(), facet));
        delete facet;
    }

    edgesList.push_back(edge_halfs[0]);
    edgesList.push_back(edge_halfs[1]);

    edgesList.erase(std::find(edgesList.begin(), edgesList.end(), this));
    delete this;
}




bool Edge3::contains(const Vertex3* vert) const
{
    if (verts[0] == vert ||
        verts[1] == vert)
        return true;

    return false;
}




bool Edge3::belongsToShell()
{
    if ((verts[0]->belongsToShellVertex ||
         verts[0]->belongsToShellEdge   ||
         verts[0]->belongsToShellFacet) &&
        (verts[1]->belongsToShellVertex ||
         verts[1]->belongsToShellEdge   ||
         verts[1]->belongsToShellFacet))
        return true;

    return false;
}




bool Edge3::needToFlip(const std::list<Facet3*>& facetsList)
{
    Vertex3* around_nodes[2];
    Facet3* around_facets[2];
    find2AdjFacets(facetsList, around_facets[0], around_facets[1]);
    around_nodes[0] = around_facets[0]->findVertNotIncludedInEdge(this);
    around_nodes[1] = around_facets[1]->findVertNotIncludedInEdge(this);

    double alpha = acos(tva::Vec3::cos(
        *verts[0] - *around_nodes[0], 
        *verts[1] - *around_nodes[0]));
    double beta = acos(tva::Vec3::cos(
        *verts[0] - *around_nodes[1], 
        *verts[1] - *around_nodes[1]));

    if (alpha + beta > PI)
        return true;

    return false;
}




Edge3::Edge3()
{
    verts[0] = nullptr;
    verts[1] = nullptr;
}


Edge3::Edge3(const Vertex3* vert0, const Vertex3* vert1)
{
    verts[0] = const_cast<Vertex3*>(vert0);
    verts[1] = const_cast<Vertex3*>(vert1);
}


Edge3::~Edge3() {}




void FrontEdge3::refreshAngleData()
{
    m_needAngleExCosProcessing = true;
    m_needComplexityProcessing = true;
}

double FrontEdge3::getComplexity()
{
    if (!m_needComplexityProcessing)
        return m_complexity;
    
    return computeComplexity();
}

double FrontEdge3::getAngleExCos()
{
    if (!m_needAngleExCosProcessing)
        return m_exCos;
    
    return computeAngleExCos();
}

double FrontEdge3::computeComplexity()
{
    m_needComplexityProcessing = false;
    return m_complexity = m_relatedCrys->getPreferredLength() / edge->magnitude() + K_ALPHA * PI / computeAngle();
}


double FrontEdge3::computeAngleExCos()
{
    auto adj_facets = getAdjFacets();

    double normals_cos = tva::Vec3::dotProduct(std::get<0>(adj_facets)->getNormal(), std::get<1>(adj_facets)->getNormal());

    m_needAngleExCosProcessing = false;
    return m_exCos = tva::spatialalgs::cpaTime(
            std::get<0>(adj_facets)->computeCenter(), std::get<0>(adj_facets)->getNormal(),
            std::get<1>(adj_facets)->computeCenter(), std::get<1>(adj_facets)->getNormal()) < 1e-6 ?
        -2.0 + normals_cos :
        -normals_cos;
}

double FrontEdge3::computeAngle()
{
    return getAngleExCos() < -1.0 ?
        acos(m_exCos + 2.0) + PI :
        acos(m_exCos);
}




std::pair<Vertex3*, Vertex3*> FrontEdge3::findOppVerts()
{
    auto adj_facets = getAdjFacets();

    return { std::get<0>(adj_facets)->facet->findVertNotIncludedInEdge(edge),
             std::get<1>(adj_facets)->facet->findVertNotIncludedInEdge(edge) };
}


FrontEdge3* FrontEdge3::findOppEdge()
{
    auto opp_verts = findOppVerts();

    for (auto& f_edge : m_relatedCrys->getFrontEdges())
    {
        if (((f_edge->edge->verts[0] == std::get<0>(opp_verts) &&
              f_edge->edge->verts[1] == std::get<1>(opp_verts)) ||
             (f_edge->edge->verts[1] == std::get<0>(opp_verts) &&
              f_edge->edge->verts[0] == std::get<1>(opp_verts))))
            return f_edge;
    }

    return nullptr;
}




std::pair<FrontFacet3*, FrontFacet3*> FrontEdge3::getAdjFacets()
{
    if (!isAdjFacetsFull())
        fillAdjFacets();

    return m_adjFacets;
}




bool FrontEdge3::addAdjFacetInPair(const FrontFacet3* fFacet)
{
    if (!std::get<0>(m_adjFacets))
        std::get<0>(m_adjFacets)  = const_cast<FrontFacet3*>(fFacet);
    else if (!std::get<1>(m_adjFacets))
        std::get<1>(m_adjFacets) = const_cast<FrontFacet3*>(fFacet);
    else
        return false;

    return true;
}


bool FrontEdge3::removeAdjFacetFromPair(const FrontFacet3* fFacet)
{
    if (std::get<0>(m_adjFacets) == const_cast<FrontFacet3*>(fFacet))
        std::get<0>(m_adjFacets) = nullptr;
    else if (std::get<1>(m_adjFacets) == const_cast<FrontFacet3*>(fFacet))
        std::get<1>(m_adjFacets) = nullptr;
    else
        return false;

    return true;
}




bool FrontEdge3::adjFacetsPairContains(const FrontFacet3* fFacet) const
{
    return (std::get<0>(m_adjFacets) == const_cast<FrontFacet3*>(fFacet)) || (std::get<1>(m_adjFacets) == const_cast<FrontFacet3*>(fFacet));
}

FrontEdge3::FrontEdge3(const Crystallite3* relatedCrys, const Edge3* edge)
    : edge(const_cast<Edge3*>(edge)), m_relatedCrys(const_cast<Crystallite3*>(relatedCrys)) {}


FrontEdge3::~FrontEdge3() {}




bool FrontEdge3::isAdjFacetsFull()
{
    return std::get<0>(m_adjFacets) && std::get<1>(m_adjFacets);
}

std::pair<FrontFacet3*, FrontFacet3*> FrontEdge3::fillAdjFacets()
{
    for (auto& f_facet : m_relatedCrys->getFrontFacets())
    {
        if (f_facet->facet->contains(edge))
        {
            if (adjFacetsPairContains(f_facet))
                continue;

            addAdjFacetInPair(f_facet);
            if (isAdjFacetsFull())
                break;
        }
    }

    return m_adjFacets;
}
