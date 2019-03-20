// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/front/surface/front-surface-face.h"
#include <stdexcept>


using FrSuFace = pmg::front::surface::Face;
using FrSuEdge = pmg::front::surface::Edge;




Vec FrSuFace::computeNormal()
{
    Vec center = computeCenter();
    Vec third_pos = face->findVertNot(face->edges[0])->pos();
    Vec loc_normal = Vec::cross(
        face->edges[0]->verts[0]->pos() - third_pos,
        face->edges[0]->verts[1]->pos() - third_pos).normalize();

    // I don't know why it led to better(correct) result.
    Vec test_normal_correct_intersect = loc_normal + Vec(static_cast<real_t>(2.1632737147),
                                                         static_cast<real_t>(1.488313178),
                                                         static_cast<real_t>(-0.71123534278))
                                                     * static_cast<real_t>(1e-3);

    int intersects_num = 0;
    for (auto& fface : m_relatedPolyhedron->frontFaces())
    {
        if (fface == this)
            continue;

        if (fface->face->intersectsBy(center, test_normal_correct_intersect))
            intersects_num++;
    }

    if (loc_normal.sqrMagnitude() < static_cast<real_t>(1e-6))
        throw std::exception();

    return normal = intersects_num % 2 == 1 ? loc_normal : -loc_normal;
}




Vec FrSuFace::computeCenter()
{
    return face->computeCenter();
}


real_t FrSuFace::computeQuality()
{
    return face->computeQuality();
}




FrSuEdge* FrSuFace::findFEdge(const pmg::Edge* edge) const
{
    for (auto& fedge : fEdges)
        if (fedge->edge == edge)
            return fedge;

    return nullptr;
}


FrSuEdge* FrSuFace::findFEdge(const pmg::Vert* v0, const pmg::Vert* v1) const
{
    for (auto& fedge : fEdges)
        if ((fedge->edge->verts[0] == v0 && fedge->edge->verts[1] == v1) ||
            (fedge->edge->verts[0] == v1 && fedge->edge->verts[1] == v0))
            return fedge;

    return nullptr;
}


FrSuEdge* FrSuFace::findFEdgeNot(const pmg::Vert* vert) const
{
    for (auto& fedge : fEdges)
        if (!fedge->edge->contains(vert))
            return fedge;

    return nullptr;
}




void FrSuFace::addFEdge(const FrSuEdge* fEdge)
{
    if (!fEdges[0])
        fEdges[0] = const_cast<FrSuEdge*>(fEdge);
    else if (!fEdges[1])
        fEdges[1] = const_cast<FrSuEdge*>(fEdge);
    else if (!fEdges[2])
        fEdges[2] = const_cast<FrSuEdge*>(fEdge);
    else
        throw std::logic_error("pmg::front::surface::Face::addFEdge can't add fEdge when fEdges already full.");
}


void FrSuFace::removeFEdge(const FrSuEdge* fEdge)
{
    for (auto& fedge : fEdges)
    {
        if (fedge == fEdge)
        {
            fedge = nullptr;
            return;
        }
    }

    throw std::logic_error("pmg::front::surface::Face::removeFEdge can't remove fEdge because there is no one.");
}


bool FrSuFace::isFEdgesFull() const
{
    if (fEdges[0] && fEdges[1] && fEdges[2])
        return true;

    return false;
}




bool FrSuFace::contains(const FrSuEdge* fEdge) const
{
    if (fEdges[0] == fEdge ||
        fEdges[1] == fEdge ||
        fEdges[2] == fEdge)
        return true;

    return false;
}




FrSuFace::Face(const Polyhedron* relatedPolyhedron, const pmg::Face* face)
    : face(const_cast<pmg::Face*>(face)), m_relatedPolyhedron(const_cast<Polyhedron*>(relatedPolyhedron)) {}
