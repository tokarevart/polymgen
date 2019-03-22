// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/face.h"
#include "helpers/spatial-algs/spatial-algs.h"

#define ONE_3 static_cast<real_t>(0.3333333333333333)


vec3 pmg::Face::computeCenter() const
{
    return ONE_3 * (
        edges[0]->verts[0]->pos() +
        edges[0]->verts[1]->pos() +
        findVertNot(edges[0])->pos());
}


real_t pmg::Face::computeQuality() const
{
    return findShortestEdge()->sqrMagnitude() / findLongestEdge()->sqrMagnitude();
}


real_t pmg::Face::computeArea() const
{
    vec3 vec0 = edges[0]->verts[1]->pos() - edges[0]->verts[0]->pos();
    vec3 vec1 = edges[1]->verts[1]->pos() - edges[1]->verts[0]->pos();
    return static_cast<real_t>(0.5) * vec3::cross(vec0, vec1).magnitude();
}




pmg::Edge* pmg::Face::intersectAlongEdge(const pmg::Face* face0, const pmg::Face* face1)
{
    int inters = 0;
    pmg::Edge* res = nullptr;
    for (auto& edge0 : face0->edges)
        for (auto& edge1 : face1->edges)
            if (edge0 == edge1)
            {
                inters++;
                res = edge0;
                break;
            }

    return inters == 1 ? res : nullptr;
}




bool pmg::Face::intersectsBy(const vec3& origin, const vec3& dir) const
{
    return spatalgs::doesRayIntersectTriangle(
        origin, dir,
        edges[0]->verts[0]->pos(),
        edges[0]->verts[1]->pos(),
        findVertNot(edges[0])->pos());
}


pmg::Vert* pmg::Face::findVertNot(const pmg::Edge* edge) const
{
    for (auto& face_edge : edges)
    {
        if (edge != face_edge)
        {
            if (!edge->contains(face_edge->verts[0]))
                return face_edge->verts[0];
            else
                return face_edge->verts[1];
        }
    }

    return nullptr;
}


pmg::Edge* pmg::Face::findEdgeNot(const pmg::Vert* vert) const
{
    for (auto& edge : edges)
    {
        if (!edge->contains(vert))
            return edge;
    }

    return nullptr;
}


pmg::Edge* pmg::Face::findEdge(const pmg::Vert* vert0, const pmg::Vert* vert1) const
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


pmg::Edge* pmg::Face::findShortestEdge() const
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


pmg::Edge* pmg::Face::findLongestEdge() const
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




bool pmg::Face::contains(const pmg::Edge* edge) const
{
    for (auto& edge0 : edges)
        if (edge0 == edge)
            return true;

    return false;
}


bool pmg::Face::contains(const pmg::Vert* vert) const
{
    for (auto& edge : edges)
        if (edge->contains(vert))
            return true;

    return false;
}




pmg::Face::Face(const pmg::Edge* edge0, const pmg::Edge* edge1, const pmg::Edge* edge2)
{
    edges[0] = const_cast<pmg::Edge*>(edge0);
    edges[1] = const_cast<pmg::Edge*>(edge1);
    edges[2] = const_cast<pmg::Edge*>(edge2);
}


pmg::Face::Face(const pmg::Vert* vert0, const pmg::Vert* vert1, const pmg::Vert* vert2)
{
    edges[0] = new pmg::Edge(vert0, vert1);
    edges[1] = new pmg::Edge(vert1, vert2);
    edges[2] = new pmg::Edge(vert2, vert0);
}
