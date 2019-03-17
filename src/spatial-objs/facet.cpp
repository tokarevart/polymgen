// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/facet.h"
#include "helpers/spatial-algs/spatial-algs.h"



tva::Vec pmg::Facet::computeCenter() const
{
    return 0.3333333333333333 * (
        edges[0]->verts[0]->pos() +
        edges[0]->verts[1]->pos() +
        findVertNot(edges[0])->pos());
}


double pmg::Facet::computeQuality() const
{
    return findShortestEdge()->sqrMagnitude() / findLongestEdge()->sqrMagnitude();
}


double pmg::Facet::computeArea() const
{
    Vec vec0 = edges[0]->verts[1]->pos() - edges[0]->verts[0]->pos();
    Vec vec1 = edges[1]->verts[1]->pos() - edges[1]->verts[0]->pos();
    return 0.5 * Vec::cross(vec0, vec1).magnitude();
}




pmg::Edge* pmg::Facet::intersectAlongEdge(const pmg::Facet* facet0, const pmg::Facet* facet1)
{
    int inters = 0;
    pmg::Edge* res = nullptr;
    for (auto& edge0 : facet0->edges)
        for (auto& edge1 : facet1->edges)
            if (edge0 == edge1)
            {
                inters++;
                res = edge0;
                break;
            }

    return inters == 1 ? res : nullptr;
}




bool pmg::Facet::intersectsBy(const tva::Point& origin, const tva::Vec& dir) const
{
    return tva::spatalgs::doesRayIntersectTriangle(
        origin, dir,
        edges[0]->verts[0]->pos(),
        edges[0]->verts[1]->pos(),
        findVertNot(edges[0])->pos());
}


pmg::Vertex* pmg::Facet::findVertNot(const pmg::Edge* edge) const
{
    for (auto& facet_edge : edges)
    {
        if (edge != facet_edge)
        {
            if (!edge->contains(facet_edge->verts[0]))
                return facet_edge->verts[0];
            else
                return facet_edge->verts[1];
        }
    }

    return nullptr;
}


pmg::Edge* pmg::Facet::findEdgeNot(const pmg::Vertex* vert) const
{
    for (auto& edge : edges)
    {
        if (!edge->contains(vert))
            return edge;
    }

    return nullptr;
}


pmg::Edge* pmg::Facet::findEdge(const pmg::Vertex* vert0, const pmg::Vertex* vert1) const
{
    for (auto& edge : edges)
    {
        if ((edge->verts[0] == vert0 &&
             edge->verts[1] == vert1) ||
            (edge->verts[0] == vert1 &&
             edge->verts[1] == vert0))
            return edge;
    }
    
    return nullptr;
}


pmg::Edge* pmg::Facet::findShortestEdge() const
{
    if (edges[0]->sqrMagnitude() < edges[1]->sqrMagnitude())
    {
        if (edges[0]->sqrMagnitude() < edges[2]->sqrMagnitude())
            return edges[0];
        else
            return edges[2];
    }
    else
    {
        if (edges[1]->sqrMagnitude() < edges[2]->sqrMagnitude())
            return edges[1];
        else
            return edges[2];
    }
}


pmg::Edge* pmg::Facet::findLongestEdge() const
{
    if (edges[0]->sqrMagnitude() > edges[1]->sqrMagnitude())
    {
        if (edges[0]->sqrMagnitude() > edges[2]->sqrMagnitude())
            return edges[0];
        else
            return edges[2];
    }
    else
    {
        if (edges[1]->sqrMagnitude() > edges[2]->sqrMagnitude())
            return edges[1];
        else
            return edges[2];
    }
}




bool pmg::Facet::contains(const pmg::Edge* edge) const
{
    for (auto& edge0 : edges)
        if (edge0 == edge)
            return true;

    return false;
}


bool pmg::Facet::contains(const pmg::Vertex* vert) const
{
    for (auto& edge : edges)
        if (edge->contains(vert))
            return true;

    return false;
}




pmg::Facet::Facet(const pmg::Edge* edge0, const pmg::Edge* edge1, const pmg::Edge* edge2)
{
    edges[0] = const_cast<pmg::Edge*>(edge0);
    edges[1] = const_cast<pmg::Edge*>(edge1);
    edges[2] = const_cast<pmg::Edge*>(edge2);
}


pmg::Facet::Facet(const pmg::Vertex* vert0, const pmg::Vertex* vert1, const pmg::Vertex* vert2)
{
    edges[0] = new pmg::Edge(vert0, vert1);
    edges[1] = new pmg::Edge(vert1, vert2);
    edges[2] = new pmg::Edge(vert2, vert0);
}
