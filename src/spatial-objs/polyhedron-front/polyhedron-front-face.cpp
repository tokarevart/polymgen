// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/polyhedron-front/polyhedron-front-face.h"
#include <stdexcept>
#include "helpers/spatial-algs/spatial-algs.h"


using namespace pmg;




vec3 front::Face::computeNormal()
{
    vec3 center = computeCenter();
    vec3 third_pos = face->findVertNot(face->edges[0])->pos();
    vec3 loc_normal = vec3::cross(
        face->edges[0]->verts[0]->pos() - third_pos,
        face->edges[0]->verts[1]->pos() - third_pos).normalize();

    // I don't know why it led to better(correct) result.
    vec3 test_normal_correct_intersect = loc_normal + vec3(static_cast<real_t>(2.1632737147),
                                                         static_cast<real_t>(1.488313178),
                                                         static_cast<real_t>(-0.71123534278))
                                                     * static_cast<real_t>(1e-3);

    size_t intersects_num = 0;
    for (auto& fface : m_relatedPolyhedron->frontFaces())
    {
        if (fface == this)
            continue;

        if (spatalgs::doesRayIntersectTriangle(
                center, test_normal_correct_intersect,
                fface->face->edges[0]->verts[0]->pos(),
                fface->face->edges[0]->verts[1]->pos(),
                fface->face->findVertNot(fface->face->edges[0])->pos()))
            intersects_num++;
    }

    if (loc_normal.sqrMagnitude() < static_cast<real_t>(1e-6))
        throw std::exception();

    return normal = intersects_num % 2 == 1 ? loc_normal : -loc_normal;
}




vec3 front::Face::computeCenter()
{
    return face->computeCenter();
}


real_t front::Face::computeQuality()
{
    return face->computeQuality();
}




front::Edge* front::Face::findFEdge(const pmg::Edge* edge) const
{
    for (auto& fedge : fEdges)
        if (fedge->edge == edge)
            return fedge;

    return nullptr;
}


front::Edge* front::Face::findFEdge(const pmg::Vert* v0, const pmg::Vert* v1) const
{
    for (auto& fedge : fEdges)
        if ((fedge->edge->verts[0] == v0 && fedge->edge->verts[1] == v1) ||
            (fedge->edge->verts[0] == v1 && fedge->edge->verts[1] == v0))
            return fedge;

    return nullptr;
}


front::Edge* front::Face::findFEdgeNot(const pmg::Vert* vert) const
{
    for (auto& fedge : fEdges)
        if (!fedge->edge->contains(vert))
            return fedge;

    return nullptr;
}




void front::Face::addFEdge(const front::Edge* fEdge)
{
    if (!fEdges[0])
        fEdges[0] = const_cast<front::Edge*>(fEdge);
    else if (!fEdges[1])
        fEdges[1] = const_cast<front::Edge*>(fEdge);
    else if (!fEdges[2])
        fEdges[2] = const_cast<front::Edge*>(fEdge);
    else
        throw std::logic_error("pmg::front::Face::addFEdge can't add fEdge when fEdges already full.");
}


void front::Face::removeFEdge(const front::Edge* fEdge)
{
    for (auto& fedge : fEdges)
    {
        if (fedge == fEdge)
        {
            fedge = nullptr;
            return;
        }
    }

    throw std::logic_error("pmg::front::Face::removeFEdge can't remove fEdge because there is no one.");
}


bool front::Face::isFEdgesFull() const
{
    if (fEdges[0] && fEdges[1] && fEdges[2])
        return true;

    return false;
}




front::Face::Face(const Polyhedron* relatedPolyhedron, const pmg::Face* face)
    : face(const_cast<pmg::Face*>(face)), m_relatedPolyhedron(const_cast<Polyhedron*>(relatedPolyhedron)) {}
