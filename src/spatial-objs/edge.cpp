// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#include "spatial-objs/edge.h"
#include <algorithm>
#include <iostream>
#include "helpers/spatial-algs/spatial-algs.h"

using pair_ff = std::pair<pmg::Facet*, pmg::Facet*>;
using tva::Vec;


#define EPS 1e-10
#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0[0], corner1[0], point[0]) && \
         BETWEEN(corner0[1], corner1[1], point[1]))


#define K_ALPHA 2.5




double pmg::Edge::magnitude() const
{
    return (*verts[1] - *verts[0]).magnitude();
}


double pmg::Edge::sqrMagnitude() const
{
    return (*verts[1] - *verts[0]).sqrMagnitude();
}




pmg::Vertex* pmg::Edge::findNot(const pmg::Edge* edge) const
{
    if (verts[0] != edge->verts[0] && verts[0] != edge->verts[1])
        return verts[0];

    if (verts[1] != edge->verts[0] && verts[1] != edge->verts[1])
        return verts[1];

    return nullptr;
}


pmg::Vertex* pmg::Edge::findNot(const pmg::Vertex* vert) const
{
    return verts[0] == vert ? verts[1] : verts[0];
}




void pmg::Edge::flip(std::list<pmg::Edge*>& edgesList, std::list<pmg::Facet*>& facetsList)
{
    pmg::Vertex* opp_nodes[2];
    auto around_facets = find2AdjFacets(facetsList);
    opp_nodes[0] = std::get<0>(around_facets)->findVertNot(this);
    opp_nodes[1] = std::get<1>(around_facets)->findVertNot(this);

    auto old_edge   = std::find(edgesList.begin(), edgesList.end(), this);
    auto old_facet0 = std::find(facetsList.begin(), facetsList.end(), std::get<0>(around_facets));
    auto old_facet1 = std::find(facetsList.begin(), facetsList.end(), std::get<1>(around_facets));

    pmg::Edge* new_edge = new pmg::Edge(opp_nodes[0], opp_nodes[1]);
    pmg::Facet* new_facet0 = new pmg::Facet(
            std::get<0>(around_facets)->findEdge(opp_nodes[0], verts[0]),
            std::get<1>(around_facets)->findEdge(opp_nodes[1], verts[0]),
            new_edge);
    pmg::Facet* new_facet1 = new pmg::Facet(
            std::get<0>(around_facets)->findEdge(opp_nodes[0], verts[1]),
            std::get<1>(around_facets)->findEdge(opp_nodes[1], verts[1]),
            new_edge);

    delete *old_edge;
    delete *old_facet0;
    delete *old_facet1;

    *old_edge = new_edge;
    *old_facet0 = new_facet0;
    *old_facet1 = new_facet1;
}


bool pmg::Edge::flipIfNeeded(std::list<pmg::Edge*>& edgesList, std::list<pmg::Facet*>& facetsList)
{
    pmg::Vertex* opp_nodes[2];
    auto around_facets = find2AdjFacets(facetsList);
    opp_nodes[0] = std::get<0>(around_facets)->findVertNot(this);
    opp_nodes[1] = std::get<1>(around_facets)->findVertNot(this);

    double alpha = acos(Vec::cos(*verts[0] - *opp_nodes[0], *verts[1] - *opp_nodes[0]));
    double beta  = acos(Vec::cos(*verts[0] - *opp_nodes[1], *verts[1] - *opp_nodes[1]));

    if (alpha + beta <= M_PI)
        return false;


    auto old_edge   = std::find(edgesList.begin(), edgesList.end(), this);
    auto old_facet0 = std::find(facetsList.begin(), facetsList.end(), std::get<0>(around_facets));
    auto old_facet1 = std::find(facetsList.begin(), facetsList.end(), std::get<1>(around_facets));

    pmg::Edge* new_edge = new pmg::Edge(opp_nodes[0], opp_nodes[1]);
    pmg::Facet* new_facet0 = new pmg::Facet(
            std::get<0>(around_facets)->findEdge(opp_nodes[0], verts[0]),
            std::get<1>(around_facets)->findEdge(opp_nodes[1], verts[0]),
            new_edge);
    pmg::Facet* new_facet1 = new pmg::Facet(
            std::get<0>(around_facets)->findEdge(opp_nodes[0], verts[1]),
            std::get<1>(around_facets)->findEdge(opp_nodes[1], verts[1]),
            new_edge);

    delete *old_edge;
    delete *old_facet0;
    delete *old_facet1;

    *old_edge = new_edge;
    *old_facet0 = new_facet0;
    *old_facet1 = new_facet1;

    return true;
}


void pmg::Edge::findAdjFacets(const std::list<pmg::Facet*>& facetsList, std::list<pmg::Facet*>& adjFacets) const
{
    for (auto& facet : facetsList)
    {
        if (facet->contains(this))
            adjFacets.push_back(facet);
    }
}


pair_ff pmg::Edge::find2AdjFacets(const std::list<pmg::Facet*>& facetsList) const
{
    pair_ff res;
    bool not_found_yet = true;
    for (auto& facet : facetsList)
    {
        if (facet->contains(this))
        {
            if (not_found_yet)
            {
                res.first = facet;
                not_found_yet = false;
            }
            else
            {
                res.second = facet;
                return res;
            }
        }
    }

    throw std::logic_error("pmg::Edge::find2AdjFacets didn't find 2 adjacent facets.");
}




bool pmg::Edge::contains(const pmg::Vertex* vert) const
{
    if (verts[0] == vert ||
        verts[1] == vert)
        return true;

    return false;
}




bool pmg::Edge::belongsToShell()
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




bool pmg::Edge::needToFlip(const std::list<pmg::Facet*>& facetsList)
{
    pmg::Vertex* opp_nodes[2];
    auto around_facets = find2AdjFacets(facetsList);
    opp_nodes[0] = std::get<0>(around_facets)->findVertNot(this);
    opp_nodes[1] = std::get<1>(around_facets)->findVertNot(this);

    double alpha = acos(Vec::cos(*verts[0] - *opp_nodes[0], *verts[1] - *opp_nodes[0]));
    double beta  = acos(Vec::cos(*verts[0] - *opp_nodes[1], *verts[1] - *opp_nodes[1]));

    if (alpha + beta <= M_PI)
        return true;

    return false;
}




pmg::Edge::Edge(const pmg::Vertex* vert0, const pmg::Vertex* vert1)
{
    verts[0] = const_cast<pmg::Vertex*>(vert0);
    verts[1] = const_cast<pmg::Vertex*>(vert1);
}
