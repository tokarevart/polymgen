// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/edge.h"
#include <algorithm>
#include <iostream>
#include "helpers/spatial-algs/spatial-algs.h"

using pair_ff = std::pair<pmg::Face*, pmg::Face*>;


#define EPS static_cast<real_t>(1e-10)
#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0[0], corner1[0], point[0]) && \
         BETWEEN(corner0[1], corner1[1], point[1]))

#define PI static_cast<real_t>(M_PI)

#define K_ALPHA static_cast<real_t>(2.5)


real_t pmg::Edge::magnitude() const
{
    return (*verts[1] - *verts[0]).magnitude();
}


real_t pmg::Edge::sqrMagnitude() const
{
    return (*verts[1] - *verts[0]).sqrMagnitude();
}


pmg::Vert* pmg::Edge::findNot(const pmg::Edge* edge) const
{
    if (verts[0] != edge->verts[0] && verts[0] != edge->verts[1])
        return verts[0];

    if (verts[1] != edge->verts[0] && verts[1] != edge->verts[1])
        return verts[1];

    return nullptr;
}


pmg::Vert* pmg::Edge::findNot(const pmg::Vert* vert) const
{
    return verts[0] == vert ? verts[1] : verts[0];
}


void pmg::Edge::flip(std::list<pmg::Edge*>& edgesList, std::list<pmg::Face*>& facesList)
{
    std::array<pmg::Vert*, 2> opp_nodes;
    auto around_faces = find2AdjFaces(facesList);
    opp_nodes[0] = std::get<0>(around_faces)->findVertNot(this);
    opp_nodes[1] = std::get<1>(around_faces)->findVertNot(this);

    auto old_edge  = std::find(edgesList.begin(), edgesList.end(), this);
    auto old_face0 = std::find(facesList.begin(), facesList.end(), std::get<0>(around_faces));
    auto old_face1 = std::find(facesList.begin(), facesList.end(), std::get<1>(around_faces));

    pmg::Edge* new_edge = new pmg::Edge(opp_nodes[0], opp_nodes[1]);
    pmg::Face* new_face0 = new pmg::Face(
            std::get<0>(around_faces)->findEdge(opp_nodes[0], verts[0]),
            std::get<1>(around_faces)->findEdge(opp_nodes[1], verts[0]),
            new_edge);
    pmg::Face* new_face1 = new pmg::Face(
            std::get<0>(around_faces)->findEdge(opp_nodes[0], verts[1]),
            std::get<1>(around_faces)->findEdge(opp_nodes[1], verts[1]),
            new_edge);

    delete *old_edge;
    delete *old_face0;
    delete *old_face1;

    *old_edge  = new_edge;
    *old_face0 = new_face0;
    *old_face1 = new_face1;
}


bool pmg::Edge::flipIfNeeded(std::list<pmg::Edge*>& edgesList, std::list<pmg::Face*>& facesList)
{
    std::array<pmg::Vert*, 2> opp_nodes;
    auto around_faces = find2AdjFaces(facesList);
    opp_nodes[0] = std::get<0>(around_faces)->findVertNot(this);
    opp_nodes[1] = std::get<1>(around_faces)->findVertNot(this);

    real_t alpha = std::acos(vec3::cos(*verts[0] - *opp_nodes[0], *verts[1] - *opp_nodes[0]));
    real_t beta  = std::acos(vec3::cos(*verts[0] - *opp_nodes[1], *verts[1] - *opp_nodes[1]));

    if (alpha + beta <= PI)
        return false;


    auto old_edge  = std::find(edgesList.begin(), edgesList.end(), this);
    auto old_face0 = std::find(facesList.begin(), facesList.end(), std::get<0>(around_faces));
    auto old_face1 = std::find(facesList.begin(), facesList.end(), std::get<1>(around_faces));

    pmg::Edge* new_edge = new pmg::Edge(opp_nodes[0], opp_nodes[1]);
    pmg::Face* new_face0 = new pmg::Face(
            std::get<0>(around_faces)->findEdge(opp_nodes[0], verts[0]),
            std::get<1>(around_faces)->findEdge(opp_nodes[1], verts[0]),
            new_edge);
    pmg::Face* new_face1 = new pmg::Face(
            std::get<0>(around_faces)->findEdge(opp_nodes[0], verts[1]),
            std::get<1>(around_faces)->findEdge(opp_nodes[1], verts[1]),
            new_edge);

    delete *old_edge;
    delete *old_face0;
    delete *old_face1;

    *old_edge  = new_edge;
    *old_face0 = new_face0;
    *old_face1 = new_face1;

    return true;
}


std::list<pmg::Face*> pmg::Edge::findAdjFaces(const std::list<pmg::Face*>& facesList) const
{
    std::list<pmg::Face*> res;
    for (auto& face : facesList)
        if (face->contains(this))
            res.push_back(face);

    return res;
}


pair_ff pmg::Edge::find2AdjFaces(const std::list<pmg::Face*>& facesList) const
{
    pair_ff res;
    bool not_found_yet = true;
    for (auto& face : facesList)
    {
        if (face->contains(this))
        {
            if (not_found_yet)
            {
                res.first = face;
                not_found_yet = false;
            }
            else
            {
                res.second = face;
                return res;
            }
        }
    }

    throw std::logic_error("pmg::Edge::find2AdjFaces didn't find 2 adjacent Faces.");
}


// TODO: replace its usage with pmg::relations content
bool pmg::Edge::contains(const pmg::Vert* vert) const
{
    if (verts[0] == vert ||
        verts[1] == vert)
        return true;

    return false;
}


// TODO: replace its usage with pmg::relations content
bool pmg::Edge::belongsToShell()
{
    if ((verts[0]->belongsToSVert ||
         verts[0]->belongsToSEdge   ||
         verts[0]->belongsToSFace) &&
        (verts[1]->belongsToSVert ||
         verts[1]->belongsToSEdge   ||
         verts[1]->belongsToSFace))
        return true;

    return false;
}


bool pmg::Edge::needToFlip(const std::list<pmg::Face*>& facesList)
{
    std::array<pmg::Vert*, 2> opp_nodes;
    auto around_faces = find2AdjFaces(facesList);
    opp_nodes[0] = std::get<0>(around_faces)->findVertNot(this);
    opp_nodes[1] = std::get<1>(around_faces)->findVertNot(this);

    real_t alpha = std::acos(vec3::cos(*verts[0] - *opp_nodes[0], *verts[1] - *opp_nodes[0]));
    real_t beta  = std::acos(vec3::cos(*verts[0] - *opp_nodes[1], *verts[1] - *opp_nodes[1]));

    if (alpha + beta <= PI)
        return true;

    return false;
}




pmg::Edge::Edge(const pmg::Vert* vert0, const pmg::Vert* vert1)
{
    verts[0] = const_cast<pmg::Vert*>(vert0);
    verts[1] = const_cast<pmg::Vert*>(vert1);
}
