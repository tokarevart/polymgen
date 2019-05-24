// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "core/front/edge.h"
#include <cmath>
#include "helpers/spatial/algs.h"


#define PI static_cast<real_t>(M_PI)

#define K_ALPHA static_cast<real_t>(4.0)


using namespace pmg;
using FrSuEdge = pmg::front::Edge;
using pair_vv = std::pair<pmg::Vert*, pmg::Vert*>;
using pair_ff = std::pair<front::Face*, front::Face*>;
using spt::vec3;


void front::Edge::refreshAngleData()
{
    m_needAngleProcessing      = true;
    m_needComplexityProcessing = true;
}


real_t front::Edge::complexity()
{
    if (!m_needComplexityProcessing)
        return m_complexity;

    return computeComplexity();
}


real_t front::Edge::angle()
{
    if (!m_needAngleProcessing)
        return m_angle;

    return computeAngle();
}


real_t front::Edge::computeComplexity()
{
    m_needComplexityProcessing = false;
    return m_complexity = m_relatedPolyhedron->preferredLength() / edge->magnitude() + K_ALPHA * PI / angle();
}


real_t front::Edge::computeAngle()
{
    auto adj_faces = adjFFaces();
    real_t normals_cos = vec3::dot(std::get<0>(adj_faces)->normal, std::get<1>(adj_faces)->normal);

    m_needAngleProcessing = false;
    return m_angle = spt::algs::cpaTime(
                std::get<0>(adj_faces)->computeCenter(), std::get<0>(adj_faces)->normal,
                std::get<1>(adj_faces)->computeCenter(), std::get<1>(adj_faces)->normal) < static_cast<real_t>(1e-6) ?
           std::acos(normals_cos) + PI :
           std::acos(-normals_cos);
}




pair_vv front::Edge::oppVerts()
{
    auto adj_faces = adjFFaces();

    return { std::get<0>(adj_faces)->face->findVertNot(edge),
             std::get<1>(adj_faces)->face->findVertNot(edge) };
}


FrSuEdge* front::Edge::findOppEdge()
{
    auto opp_verts = oppVerts();
    std::vector<front::Edge*> opp_fedges;
    for (auto& f_edge : m_relatedPolyhedron->frontEdges())
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

    std::vector<pair_ff> adj_ffaces_vec;
    adj_ffaces_vec.reserve(opp_fedges.size());
    for (auto& fedge : opp_fedges)
        adj_ffaces_vec.push_back(fedge->adjFFaces());

    vec3 main_vert_pos = edge->verts[0]->pos();
    vec3 main_vert_proj = spt::algs::project(main_vert_pos, opp_verts.first->pos(), opp_verts.second->pos());
    vec3 main_vec = main_vert_pos - main_vert_proj;

    front::Edge* max_cos_fedge = nullptr;
    real_t max_cos = -1.0;
    for (std::size_t i = 0; i < adj_ffaces_vec.size(); i++)
    {
        vec3 adj_opp_pos0 = adj_ffaces_vec[i].first->face->findVertNot(opp_fedges.front()->edge)->pos();
        vec3 adj_opp_pos1 = adj_ffaces_vec[i].second->face->findVertNot(opp_fedges.front()->edge)->pos();
        vec3 adj_opp_proj0 = spt::algs::project(adj_opp_pos0, opp_verts.first->pos(), opp_verts.second->pos());
        vec3 adj_opp_proj1 = spt::algs::project(adj_opp_pos1, opp_verts.first->pos(), opp_verts.second->pos());
        vec3 adj_vec0 = adj_opp_pos0 - adj_opp_proj0;
        vec3 adj_vec1 = adj_opp_pos1 - adj_opp_proj1;
        real_t cos0 = vec3::cos(adj_vec0, main_vec);
        real_t cos1 = vec3::cos(adj_vec1, main_vec);
        real_t max_cur_cos = std::max(cos0, cos1);
        if (max_cur_cos > max_cos)
        {
            max_cos = max_cur_cos;
            max_cos_fedge = opp_fedges[i];
        }
    }

    return max_cos_fedge;
}




pair_ff front::Edge::adjFFaces()
{
    if (!adjFacesFull())
        fillAdjFFaces();

    return m_adjFFaces;
}




bool front::Edge::addAdjFFace(const front::Face* fFace)
{
    if (!std::get<0>(m_adjFFaces))
        std::get<0>(m_adjFFaces)  = const_cast<front::Face*>(fFace);
    else if (!std::get<1>(m_adjFFaces))
        std::get<1>(m_adjFFaces) = const_cast<front::Face*>(fFace);
    else
        return false;

    return true;
}


bool front::Edge::removeAdjFFace(const front::Face* fFace)
{
    if (std::get<0>(m_adjFFaces) == const_cast<front::Face*>(fFace))
        std::get<0>(m_adjFFaces) = nullptr;
    else if (std::get<1>(m_adjFFaces) == const_cast<front::Face*>(fFace))
        std::get<1>(m_adjFFaces) = nullptr;
    else
        return false;

    return true;
}


bool front::Edge::adjFFacesContains(const front::Face* fFace) const
{
    return std::get<0>(m_adjFFaces) == const_cast<front::Face*>(fFace) ||
           std::get<1>(m_adjFFaces) == const_cast<front::Face*>(fFace);
}


void front::Edge::fillAdjFFaces(const front::Face* fFace0, const front::Face* fFace1)
{
    m_adjFFaces.first  = const_cast<front::Face*>(fFace0);
    m_adjFFaces.second = const_cast<front::Face*>(fFace1);
}




front::Edge::Edge(const Polyhedron* relatedPolyhedron, const pmg::Edge* edge)
    : edge(const_cast<pmg::Edge*>(edge)), m_relatedPolyhedron(const_cast<Polyhedron*>(relatedPolyhedron)) {}




bool front::Edge::adjFacesFull()
{
    return std::get<0>(m_adjFFaces) && std::get<1>(m_adjFFaces);
}


pair_ff front::Edge::fillAdjFFaces()
{
    for (auto& fface : m_relatedPolyhedron->frontFaces())
    {
        if (fface->contains(this))
        {
            if (adjFFacesContains(fface))
                continue;

            addAdjFFace(fface);
            if (adjFacesFull())
                break;
        }
    }

    if (!adjFacesFull())
        throw std::logic_error("pmg::front::Edge::fillAdjFFaces didn't find 2 adjacent front Faces.");

    return m_adjFFaces;
}
