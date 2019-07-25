// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "edge.h"
#include <algorithm>
#include <iostream>
#include "../helpers/mathconsts.h"
#include "../helpers/spatial/algs.h"

using pair_ff = std::pair<pmg::Face*, pmg::Face*>;
using vec3 = spt::vec<3, real_t>;


real_t pmg::Edge::magnitude() const {
    return (verts[1]->pos() - verts[0]->pos()).magnitude();
}


real_t pmg::Edge::sqr_magnitude() const {
    return (verts[1]->pos() - verts[0]->pos()).sqr_magnitude();
}


pmg::Vert* pmg::Edge::vert_not(const pmg::Edge* edge) const {
    if (verts[0] != edge->verts[0] && verts[0] != edge->verts[1])
        return verts[0];

    if (verts[1] != edge->verts[0] && verts[1] != edge->verts[1])
        return verts[1];

    return nullptr;
}


pmg::Vert* pmg::Edge::vert_not(const pmg::Vert* vert) const {
    return verts[0] == vert ? verts[1] : verts[0];
}


void pmg::Edge::flip(std::list<pmg::Edge*>& edgesList, std::list<pmg::Face*>& facesList) {
    std::array<pmg::Vert*, 2> opp_nodes;
    auto around_faces = adj_2_faces(facesList);
    opp_nodes[0] = std::get<0>(around_faces)->find_vert_not(this);
    opp_nodes[1] = std::get<1>(around_faces)->find_vert_not(this);

    auto old_edge = std::find(edgesList.begin(), edgesList.end(), this);
    auto old_face0 = std::find(facesList.begin(), facesList.end(), std::get<0>(around_faces));
    auto old_face1 = std::find(facesList.begin(), facesList.end(), std::get<1>(around_faces));

    pmg::Edge* new_edge = new pmg::Edge(opp_nodes[0], opp_nodes[1]);
    pmg::Face* new_face0 = new pmg::Face(
        std::get<0>(around_faces)->find_edge(opp_nodes[0], verts[0]),
        std::get<1>(around_faces)->find_edge(opp_nodes[1], verts[0]),
        new_edge);
    pmg::Face* new_face1 = new pmg::Face(
        std::get<0>(around_faces)->find_edge(opp_nodes[0], verts[1]),
        std::get<1>(around_faces)->find_edge(opp_nodes[1], verts[1]),
        new_edge);

    delete *old_edge;
    delete *old_face0;
    delete *old_face1;

    *old_edge = new_edge;
    *old_face0 = new_face0;
    *old_face1 = new_face1;
}


bool pmg::Edge::flip_if_needed(std::list<pmg::Edge*>& edgesList, std::list<pmg::Face*>& facesList) {
    std::array<pmg::Vert*, 2> opp_nodes;
    auto around_faces = adj_2_faces(facesList);
    opp_nodes[0] = std::get<0>(around_faces)->find_vert_not(this);
    opp_nodes[1] = std::get<1>(around_faces)->find_vert_not(this);
    
    real_t alpha = std::acos(spt::cos(verts[0]->pos() - opp_nodes[0]->pos(), verts[1]->pos() - opp_nodes[0]->pos()));
    real_t beta = std::acos(spt::cos(verts[0]->pos() - opp_nodes[1]->pos(), verts[1]->pos() - opp_nodes[1]->pos()));

    if (alpha + beta <= mathconsts::PI)
        return false;


    auto old_edge = std::find(edgesList.begin(), edgesList.end(), this);
    auto old_face0 = std::find(facesList.begin(), facesList.end(), std::get<0>(around_faces));
    auto old_face1 = std::find(facesList.begin(), facesList.end(), std::get<1>(around_faces));

    pmg::Edge* new_edge = new pmg::Edge(opp_nodes[0], opp_nodes[1]);
    pmg::Face* new_face0 = new pmg::Face(
        std::get<0>(around_faces)->find_edge(opp_nodes[0], verts[0]),
        std::get<1>(around_faces)->find_edge(opp_nodes[1], verts[0]),
        new_edge);
    pmg::Face* new_face1 = new pmg::Face(
        std::get<0>(around_faces)->find_edge(opp_nodes[0], verts[1]),
        std::get<1>(around_faces)->find_edge(opp_nodes[1], verts[1]),
        new_edge);

    delete *old_edge;
    delete *old_face0;
    delete *old_face1;

    *old_edge = new_edge;
    *old_face0 = new_face0;
    *old_face1 = new_face1;

    return true;
}


std::list<pmg::Face*> pmg::Edge::adj_faces(const std::list<pmg::Face*>& facesList) const {
    std::list<pmg::Face*> res;
    for (auto& face : facesList)
        if (face->contains(this))
            res.push_back(face);

    return res;
}


pair_ff pmg::Edge::adj_2_faces(const std::list<pmg::Face*>& facesList) const {
    pair_ff res;
    bool not_found_yet = true;
    for (auto& face : facesList) {
        if (face->contains(this)) {
            if (not_found_yet) {
                res.first = face;
                not_found_yet = false;
            } else {
                res.second = face;
                return res;
            }
        }
    }

    throw std::logic_error("pmg::Edge::adj_2_faces didn't find 2 adjacent Faces.");
}


bool pmg::Edge::contains(const pmg::Vert* vert) const {
    if (verts[0] == vert ||
        verts[1] == vert)
        return true;

    return false;
}


bool pmg::Edge::belongs_to_shell() {
    if ((verts[0]->belongs_to_svert ||
         verts[0]->belongs_to_sedge ||
         verts[0]->belongs_to_sface) &&
         (verts[1]->belongs_to_svert ||
          verts[1]->belongs_to_sedge ||
          verts[1]->belongs_to_sface))
        return true;

    return false;
}


bool pmg::Edge::need_to_flip(const std::list<pmg::Face*>& facesList) {
    std::array<pmg::Vert*, 2> opp_nodes;
    auto around_faces = adj_2_faces(facesList);
    opp_nodes[0] = std::get<0>(around_faces)->find_vert_not(this);
    opp_nodes[1] = std::get<1>(around_faces)->find_vert_not(this);

    real_t alpha = std::acos(spt::cos(verts[0]->pos() - opp_nodes[0]->pos(), verts[1]->pos() - opp_nodes[0]->pos()));
    real_t beta = std::acos(spt::cos(verts[0]->pos() - opp_nodes[1]->pos(), verts[1]->pos() - opp_nodes[1]->pos()));

    if (alpha + beta <= mathconsts::PI)
        return true;

    return false;
}




pmg::Edge::Edge(const pmg::Vert* vert0, const pmg::Vert* vert1) {
    verts[0] = const_cast<pmg::Vert*>(vert0);
    verts[1] = const_cast<pmg::Vert*>(vert1);
}
