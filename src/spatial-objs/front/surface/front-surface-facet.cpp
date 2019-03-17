// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/front/surface/front-surface-facet.h"
#include <stdexcept>


using FrSuFacet = pmg::front::surface::Facet;
using FrSuEdge  = pmg::front::surface::Edge;
using tva::Vec;




Vec FrSuFacet::computeNormal()
{
    Vec center = computeCenter();
    Vec third_pos = facet->findVertNot(facet->edges[0])->pos();
    Vec loc_normal = Vec::cross(
        facet->edges[0]->verts[0]->pos() - third_pos,
        facet->edges[0]->verts[1]->pos() - third_pos).normalize();

    // I don't know why it led to better(correct) result.
    Vec test_normal_correct_intersect = loc_normal + Vec(2.1632737147, 1.488313178, -0.71123534278) * 1e-3;

    int intersects_num = 0;
    for (auto& f_facet : m_relatedCrys->frontFacets())
    {
        if (f_facet == this)
            continue;

        if (f_facet->facet->intersectsBy(center, test_normal_correct_intersect))
            intersects_num++;
    }

    if (loc_normal.sqrMagnitude() < 1e-6)
        throw std::exception();

    return normal = intersects_num % 2 == 1 ? loc_normal : -loc_normal;
}




Vec FrSuFacet::computeCenter()
{
    return facet->computeCenter();
}


double FrSuFacet::computeQuality()
{
    return facet->computeQuality();
}




FrSuEdge* FrSuFacet::findFEdge(const pmg::Edge* edge) const
{
    for (auto& fedge : fEdges)
        if (fedge->edge == edge)
            return fedge;

    return nullptr;
}


FrSuEdge* FrSuFacet::findFEdge(const pmg::Vertex* v0, const pmg::Vertex* v1) const
{
    for (auto& fedge : fEdges)
        if ((fedge->edge->verts[0] == v0 && fedge->edge->verts[1] == v1) ||
            (fedge->edge->verts[0] == v1 && fedge->edge->verts[1] == v0))
            return fedge;

    return nullptr;
}


FrSuEdge* FrSuFacet::findFEdgeNot(const pmg::Vertex* vert) const
{
    for (auto& fedge : fEdges)
        if (!fedge->edge->contains(vert))
            return fedge;

    return nullptr;
}




void FrSuFacet::addFEdge(const FrSuEdge* fEdge)
{
    if (!fEdges[0])
        fEdges[0] = const_cast<FrSuEdge*>(fEdge);
    else if (!fEdges[1])
        fEdges[1] = const_cast<FrSuEdge*>(fEdge);
    else if (!fEdges[2])
        fEdges[2] = const_cast<FrSuEdge*>(fEdge);
    else
        throw std::logic_error("pmg::front::surface::Facet::addFEdge can't add fEdge when fEdges already full.");
}


void FrSuFacet::removeFEdge(const FrSuEdge* fEdge)
{
    for (auto& fedge : fEdges)
    {
        if (fedge == fEdge)
        {
            fedge = nullptr;
            return;
        }
    }

    throw std::logic_error("pmg::front::surface::Facet::removeFEdge can't remove fEdge because there is no one.");
}


bool FrSuFacet::isFEdgesFull() const
{
    if (fEdges[0] && fEdges[1] && fEdges[2])
        return true;

    return false;
}




bool FrSuFacet::contains(const FrSuEdge* fEdge) const
{
    if (fEdges[0] == fEdge ||
        fEdges[1] == fEdge ||
        fEdges[2] == fEdge)
        return true;

    return false;
}




FrSuFacet::Facet(const Crystallite* relatedCrys, const pmg::Facet* facet)
    : facet(const_cast<pmg::Facet*>(facet)), m_relatedCrys(const_cast<Crystallite*>(relatedCrys)) {}
