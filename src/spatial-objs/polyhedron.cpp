// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "polyhedron.h"
#include <algorithm>
#include <iostream>
#include <float.h>
#include "helpers/spatial-algs/spatial-algs.h"

#include "helpers/cosd-values.h"

using namespace pmg;

using pair_ff = std::pair<front::Face*, front::Face*>;


#define DET(a, b, c, d) \
        (a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0[0], corner1[0], point[0]) && \
         BETWEEN(corner0[1], corner1[1], point[1]))

#define ALPHA_P      static_cast<real_t>(70.52877936550931)
#define DEG_1_IN_RAD static_cast<real_t>( 0.0174532925199432957)

#define ONE_3                static_cast<real_t>(0.3333333333333333)
#define SQRT3_2              static_cast<real_t>(0.8660254037844386)
#define SQRT_2_3             static_cast<real_t>(0.8164965809277260)
#define ONE_PLUS_SQRT2_SQRT3 static_cast<real_t>(1.3938468501173517)

#define NOT_TOO_CLOSE          static_cast<real_t>(2e-1)
#define FROM_VERT_COEF         static_cast<real_t>(1e-2)
#define EDGES_INTERS_DIST_COEF static_cast<real_t>(4e-3)

#define K_MAXD static_cast<real_t>(0.3)
#define K_D    static_cast<real_t>(0.4)


template <typename T>
constexpr real_t degToRad( T value )
{
    return value * DEG_1_IN_RAD;
}


//#define DEV_DEBUG




bool Polyhedron::shellContains(const Vert* vert) const
{
    if (vert->belongsToShellFace)
    {
        for (auto& sface : m_shellFaces)
            if (sface == vert->belongsToShellFace)
                return true;
    }
    else if (vert->belongsToShellEdge)
    {
        for (auto& sedge : m_shellEdges)
            if (sedge == vert->belongsToShellEdge)
                return true;
    }
    else if (vert->belongsToShellVert)
    {
        for (auto& sedge : m_shellEdges)
            if (sedge->verts[0] == vert->belongsToShellVert ||
                sedge->verts[1] == vert->belongsToShellVert)
                return true;
    }

    return false;
}




void Polyhedron::initializeFFaceFEdges(front::Face* fFace) const
{
    int n_added = 0;
    for (auto& fedge : m_frontEdges)
    {
        if (fFace->face->contains(fedge->edge))
        {
            fFace->addFEdge(fedge);
            if (++n_added == 3)
                return;
        }
    }

    throw std::logic_error("Polyhedron::initializeFFaceFEdges didn't find 3 front edges.");
}


void Polyhedron::initializeFront()
{
    for (auto& sedge : m_shellEdges)
        for (auto& edge : sedge->innerEdges())
            m_frontEdges.push_back(new front::Edge(this, edge));

    for (auto& sface : m_shellFaces)
        for (auto& edge : sface->innerEdges())
            m_frontEdges.push_back(new front::Edge(this, edge));

    for (auto& sface : m_shellFaces)
        for (auto& face : sface->innerFaces())
            m_frontFaces.push_back(new front::Face(this, face));

    for (auto& fface : m_frontFaces)
        initializeFFaceFEdges(fface);
}


void Polyhedron::computeFrontNormals()
{
    for (auto& face : m_frontFaces)
        face->computeNormal();
}




shell::Edge* Polyhedron::findShellEdge(const shell::Vert* v0, const shell::Vert* v1) const
{
    for (auto& sedge : m_shellEdges)
    {
        if ((sedge->verts[0] == v0  &&
             sedge->verts[1] == v1) ||
            (sedge->verts[1] == v0  &&
             sedge->verts[0] == v1))
            return sedge;
    }

    return nullptr;
}


front::Face* Polyhedron::findFrontFace(const Face* face) const
{
    for (auto& fface : m_frontFaces)
    {
        if (fface->face == face)
            return fface;
    }

    return nullptr;
}


std::vector<front::Edge*> Polyhedron::findFEdge(const Vert* v0, const Vert* v1) const
{
    std::vector<front::Edge*> res;
    for (auto& f_edge : m_frontEdges)
    {
        if (((f_edge->edge->verts[0] == v0 &&
              f_edge->edge->verts[1] == v1) ||
             (f_edge->edge->verts[1] == v0 &&
              f_edge->edge->verts[0] == v1)))
        {
            res.push_back(f_edge);
        }
    }

    return res;
}


std::vector<front::Edge*> Polyhedron::findFEdge(const Edge* edge) const
{
    std::vector<front::Edge*> res;
    for (auto& f_edge : m_frontEdges)
    {
        if (f_edge->edge == edge)
        {
            res.push_back(f_edge);
        }
    }

    return res;
}




front::Face* Polyhedron::addToFront(const Face* face, bool addInner)
{
    front::Face* new_fface = new front::Face(this, face);
    m_frontFaces.push_back(new_fface);
    if (addInner)
        m_innerFaces.push_back(new_fface->face);
    return new_fface;
}


front::Edge* Polyhedron::addToFront(const pmg::Edge* edge, bool addInner)
{
    front::Edge* new_f_edge = new front::Edge(this, edge);
    m_frontEdges.push_back(new_f_edge);
    if (addInner)
        m_innerEdges.push_back(new_f_edge->edge);
    return new_f_edge;
}




void Polyhedron::removeFromFront(front::Face* fFace)
{
    m_frontFaces.erase(std::find(m_frontFaces.begin(), m_frontFaces.end(), fFace));
    delete fFace;
}


void Polyhedron::removeFromFront(front::Edge* fEdge)
{
    m_frontEdges.erase(std::find(m_frontEdges.begin(), m_frontEdges.end(), fEdge));
    delete fEdge;
}




bool Polyhedron::vertInsideFrontCheck(const Vec& v) const
{
    int inters_num = 0;
    Vec dir = ONE_3 * (m_frontFaces.front()->computeCenter() - static_cast<real_t>(3.0) * v);

    for (auto& fface : m_frontFaces)
    {
        if (spatalgs::doesRayIntersectTriangle(
                v, dir,
                fface->face->edges[0]->verts[0]->pos(),
                fface->face->edges[0]->verts[1]->pos(),
                fface->face->findVertNot(fface->face->edges[0])->pos()))
        {
             inters_num++;
        }
    }

    return inters_num % 2 == 1;
}


bool Polyhedron::segmentGlobalIntersectionCheck(const Vec& v0, const Vec& v1) const
{
    for (auto& face : m_innerFaces)
    {
        if (spatalgs::doesSegmentIntersectTriangle(
                v0, v1,
                face->edges[0]->verts[0]->pos(),
                face->edges[0]->verts[1]->pos(),
                face->findVertNot(face->edges[0])->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::segmentFrontIntersectionCheck(const Vec& v0, const Vec& v1) const
{
    for (auto& fface : m_frontFaces)
    {
        if (spatalgs::doesSegmentIntersectTriangle(
                v0, v1,
                fface->face->edges[0]->verts[0]->pos(),
                fface->face->edges[0]->verts[1]->pos(),
                fface->face->findVertNot(fface->face->edges[0])->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::edgeGlobalIntersectionCheck(const Edge* edge) const
{
    Vec delta = static_cast<real_t>(8.0) * FROM_VERT_COEF * (*edge->verts[1] - *edge->verts[0]);
    Vec segment[2];
    segment[0] = edge->verts[0]->pos() + delta;
    segment[1] = edge->verts[1]->pos() - delta;

    for (auto& face : m_innerFaces)
    {
        if (!face->contains(edge) &&
            spatalgs::doesSegmentIntersectTriangle(
                segment[0], segment[1],
                face->edges[0]->verts[0]->pos(),
                face->edges[0]->verts[1]->pos(),
                face->findVertNot(face->edges[0])->pos()))
            return true;
    }

    return false;
}


bool XOR(bool b0, bool b1)
{
    return (b0 || b1) && !(b0 && b1);
}


bool Polyhedron::edgeIntersectionCheck(front::Edge* frontEdge) const
{
    auto opp_verts = frontEdge->findOppVerts();
    Vec opp_verts_poses[2];
    opp_verts_poses[0] = std::get<0>(opp_verts)->pos();
    opp_verts_poses[1] = std::get<1>(opp_verts)->pos();
    Vec opp_vert0_to_1 = opp_verts_poses[1] - opp_verts_poses[0];
    Vec delta_vec_0th_to_1st = FROM_VERT_COEF * opp_vert0_to_1;

    if (segmentFrontIntersectionCheck(
        opp_verts_poses[0] + delta_vec_0th_to_1st,
        opp_verts_poses[1] - delta_vec_0th_to_1st))
        return true;

    for (auto& f_edge : m_frontEdges)
    {
        if (f_edge->edge->contains(std::get<0>(opp_verts)) &&
            f_edge->edge->contains(std::get<1>(opp_verts)))
            return false;
    }

    for (auto& f_edge : m_frontEdges)
    {
        Vert* vert_buf;
        if (bool contains[2] { f_edge->edge->contains(std::get<0>(opp_verts)),
                               f_edge->edge->contains(std::get<1>(opp_verts)) };
            contains[0])
        {
            if (f_edge->edge->verts[0] == std::get<0>(opp_verts))
                vert_buf = f_edge->edge->verts[1];
            else
                vert_buf = f_edge->edge->verts[0];

            if (spatalgs::distancePointToSegment(vert_buf->pos(), opp_verts_poses[0], opp_verts_poses[1]) < EDGES_INTERS_DIST_COEF * m_preferredLength)
                return true;
        }
        else if (contains[1])
        {
            if (f_edge->edge->verts[0] == std::get<1>(opp_verts))
                vert_buf = f_edge->edge->verts[1];
            else
                vert_buf = f_edge->edge->verts[0];

            if (spatalgs::distancePointToSegment(vert_buf->pos(), opp_verts_poses[0], opp_verts_poses[1]) < EDGES_INTERS_DIST_COEF * m_preferredLength)
                return true;
        }
        else
        {
            if (spatalgs::segmentsDistance(
                    opp_verts_poses[0], opp_verts_poses[1],
                    f_edge->edge->verts[0]->pos(), f_edge->edge->verts[1]->pos()) < EDGES_INTERS_DIST_COEF * m_preferredLength)
                return true;
        }
    }

    return false;
}


bool Polyhedron::faceIntersectionCheck(const Vert* v0, const Vert* v1, const Vec& v2) const
{
    for (auto& f_edge : m_frontEdges)
    {
        bool contains[2];
        contains[0] = f_edge->edge->contains(v0);
        contains[1] = f_edge->edge->contains(v1);

        if (contains[0] || contains[1])
            continue;

        Vec first_to_second_delta = NOT_TOO_CLOSE * (f_edge->edge->verts[1]->pos() - f_edge->edge->verts[0]->pos());
        Vec f_edge_verts_poses[2];
        f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos() - first_to_second_delta;
        f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos() + first_to_second_delta;

        if (spatalgs::doesSegmentIntersectTriangle(
            f_edge_verts_poses[0], f_edge_verts_poses[1],
            v0->pos(), v1->pos(), v2))
            return true;
    }

    return false;
}


bool Polyhedron::faceIntersectionCheck(const Vert* v0, const Vert* v1, const Vert* v2) const
{
    Vec verts_poses[3];
    verts_poses[0] = v0->pos();
    verts_poses[1] = v1->pos();
    verts_poses[2] = v2->pos();

    for (auto& f_edge : m_frontEdges)
    {
        bool contains[3];
        contains[0] = f_edge->edge->contains(v0);
        contains[1] = f_edge->edge->contains(v1);
        contains[2] = f_edge->edge->contains(v2);

        if ((contains[0] && contains[1]) ||
            (contains[0] && contains[2]) ||
            (contains[1] && contains[2]))
            continue;

        Vec first_to_second_delta = FROM_VERT_COEF * (f_edge->edge->verts[1]->pos() - f_edge->edge->verts[0]->pos());
        Vec f_edge_verts_poses[2];
        if (contains[0])
        {
            if (f_edge->edge->verts[0] == v0)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos() - first_to_second_delta;
            }
        }
        else if (contains[1])
        {
            if (f_edge->edge->verts[0] == v1)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos() - first_to_second_delta;
            }
        }
        else if (contains[2])
        {
            if (f_edge->edge->verts[0] == v2)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos() - first_to_second_delta;
            }
        }
        else
        {
            f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos();
            f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos();
        }

        if (spatalgs::doesSegmentIntersectTriangle(
                f_edge_verts_poses[0], f_edge_verts_poses[1],
                verts_poses[0], verts_poses[1], verts_poses[2]))
            return true;
    }

    return false;
}


bool Polyhedron::facesIntersectionCheck(front::Edge* frontEdge) const
{
    auto opp_verts = frontEdge->findOppVerts();
        
    if (faceIntersectionCheck(frontEdge->edge->verts[0], std::get<0>(opp_verts), std::get<1>(opp_verts)) ||
        faceIntersectionCheck(frontEdge->edge->verts[1], std::get<0>(opp_verts), std::get<1>(opp_verts)))
        return true;

    return false;
}


bool Polyhedron::insideTetrCheck(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3, const Vec& vert) const
{
    Vec vert_to_p0 = p0 - vert;
    Vec vert_to_p1 = p1 - vert;
    Vec vert_to_p2 = p2 - vert;
    Vec vert_to_p3 = p3 - vert;

    real_t abs_mixed_prods[5];
    abs_mixed_prods[0] = std::abs(Vec::mixed(vert_to_p0, vert_to_p2, vert_to_p3));
    abs_mixed_prods[1] = std::abs(Vec::mixed(vert_to_p0, vert_to_p1, vert_to_p2));
    abs_mixed_prods[2] = std::abs(Vec::mixed(vert_to_p0, vert_to_p1, vert_to_p3));
    abs_mixed_prods[3] = std::abs(Vec::mixed(vert_to_p1, vert_to_p2, vert_to_p3));
    abs_mixed_prods[4] = std::abs(Vec::mixed(p1 - p0, p2 - p0, p3 - p0));

    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4]/* * 1.000001*/;
}


bool Polyhedron::anyVertInsidePotentialTetrCheck(front::Edge* fEdge) const
{
    auto opp_verts = fEdge->findOppVerts();

    Vec points[4];
    points[0] = std::get<0>(opp_verts)->pos();
    points[1] = std::get<1>(opp_verts)->pos();
    points[2] = fEdge->edge->verts[0]->pos();
    points[3] = fEdge->edge->verts[1]->pos();
    for (auto& vert : m_innerVerts)
    {
        if (vert != std::get<0>(opp_verts) &&
            vert != std::get<1>(opp_verts) &&
            vert != fEdge->edge->verts[0] &&
            vert != fEdge->edge->verts[1] &&
            insideTetrCheck(points[0], points[1], points[2], points[3], vert->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::frontSplitCheck(front::Edge* fEdge, front::Edge* oppFEdge) const
{
    auto opp_fedge = oppFEdge;
    if (!opp_fedge)
        opp_fedge = fEdge->findOppEdge();
    if (!opp_fedge)
        return false;

    auto opp_fedge_opp_verts = opp_fedge->findOppVerts();

    if (std::get<0>(opp_fedge_opp_verts) == fEdge->edge->verts[0] ||
        std::get<0>(opp_fedge_opp_verts) == fEdge->edge->verts[1] ||
        std::get<1>(opp_fedge_opp_verts) == fEdge->edge->verts[0] ||
        std::get<1>(opp_fedge_opp_verts) == fEdge->edge->verts[1])
        return false;

    return true;
}


bool Polyhedron::frontCollapseCheck(front::Edge* fEdge, front::Edge* oppFEdge) const
{
    auto opp_fedge = oppFEdge;
    if (!opp_fedge)
        opp_fedge = fEdge->findOppEdge();
    if (!opp_fedge)
        return false;

    auto opp_ffaces = opp_fedge->getAdjFFaces();

    if (!fEdge->edge->contains(opp_ffaces.first->face->findVertNot(opp_fedge->edge)) ||
        !fEdge->edge->contains(opp_ffaces.second->face->findVertNot(opp_fedge->edge)))
        return false;

    return true;
}


bool Polyhedron::parallelFacesCheck(front::Edge* fEdge) const
{
    auto adj_ffaces = fEdge->getAdjFFaces();

    Vert* opp_verts[2];
    opp_verts[0] = std::get<0>(adj_ffaces)->face->findVertNot(fEdge->edge);
    opp_verts[1] = std::get<1>(adj_ffaces)->face->findVertNot(fEdge->edge);

    Vec plane0[2];
    plane0[0] = *opp_verts[0] - *fEdge->edge->verts[0];
    plane0[1] = *opp_verts[1] - *fEdge->edge->verts[0];
    Vec normal0 = Vec::cross(plane0[0], plane0[1]).normalize();
    Vec plane1[2];
    plane1[0] = *opp_verts[0] - *fEdge->edge->verts[1];
    plane1[1] = *opp_verts[1] - *fEdge->edge->verts[1];
    Vec normal1 = Vec::cross(plane1[0], plane1[1]).normalize();

    for (auto& fface : m_frontFaces)
    {
        Edge* inters_reses[2];
        if ((fface != std::get<0>(adj_ffaces)) &&
            (fface != std::get<1>(adj_ffaces)))
        {
            inters_reses[0] = Face::intersectAlongEdge(fface->face, std::get<0>(adj_ffaces)->face);
            inters_reses[1] = Face::intersectAlongEdge(fface->face, std::get<1>(adj_ffaces)->face);
            if (!XOR(static_cast<bool>(inters_reses[0]), static_cast<bool>(inters_reses[1])))
                continue;

            Vert* fface_to_verts[2];
            fface_to_verts[0] = fface->face->edges[0]->verts[0];
            fface_to_verts[1] = fface->face->edges[0]->verts[1];
            Vert* fface_from_vert = fface->face->findVertNot(fface->face->edges[0]);

            Vec f_plane[2];
            f_plane[0] = *fface_to_verts[0] - *fface_from_vert;
            f_plane[1] = *fface_to_verts[1] - *fface_from_vert;
            Vec f_normal = Vec::cross(f_plane[0], f_plane[1]).normalize();

            if (std::abs(std::abs(Vec::dot(f_normal, normal0)) - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6) ||
                std::abs(std::abs(Vec::dot(f_normal, normal1)) - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6))
            {
                int i = inters_reses[0] ? 0 : 1;

                Vec border_verts[2];
                border_verts[0] = inters_reses[i]->verts[0]->pos();
                border_verts[1] = inters_reses[i]->verts[1]->pos();

                Vec main_face_3rd_vert;
                if (inters_reses[i]->contains(opp_verts[0]))
                    main_face_3rd_vert = opp_verts[1]->pos();
                else
                    main_face_3rd_vert = opp_verts[0]->pos();

                Vec cur_face_3rd_vert = fface->face->findVertNot(inters_reses[i])->pos();

                Vec main_face_cross = Vec::cross(main_face_3rd_vert - border_verts[0], main_face_3rd_vert - border_verts[1]);
                Vec cur_face_cross = Vec::cross(cur_face_3rd_vert - border_verts[0], cur_face_3rd_vert - border_verts[1]);
                if (Vec::dot(main_face_cross, cur_face_cross) > static_cast<real_t>(0.0))
                    return true;
            }
        }
    }

    return false;
}


bool Polyhedron::doesFrontIntersectSphere(const Vec& center, real_t radius) const
{
    for (auto& fface : m_frontFaces)
    {
        Vec triangle[3];
        triangle[0] = fface->face->edges[0]->verts[0]->pos();
        triangle[1] = fface->face->edges[0]->verts[1]->pos();
        triangle[2] = fface->face->findVertNot(fface->face->edges[0])->pos();
        if (spatalgs::doesTriangleIntersectSphere(triangle[0], triangle[1], triangle[2], center, radius))
            return true;
    }

    return false;
}




std::pair<real_t, real_t> Polyhedron::computeMinMaxEdgesLengths(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3)
{
    auto min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    min_max.first  = sqrtReal(min_max.first);
    min_max.second = sqrtReal(min_max.second);
    return min_max;
}


std::pair<real_t, real_t> Polyhedron::computeMinMaxEdgesSqrLengths(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3)
{
    real_t sqr_magns[6];
    sqr_magns[0] = (p1 - p0).sqrMagnitude();
    sqr_magns[1] = (p2 - p0).sqrMagnitude();
    sqr_magns[2] = (p3 - p0).sqrMagnitude();
    sqr_magns[3] = (p2 - p1).sqrMagnitude();
    sqr_magns[4] = (p3 - p1).sqrMagnitude();
    sqr_magns[5] = (p3 - p2).sqrMagnitude();
    return std::minmax({ sqr_magns[0], sqr_magns[1], sqr_magns[2], sqr_magns[3], sqr_magns[4], sqr_magns[5] });
}


real_t Polyhedron::computeTetrSimpleQuality(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return sqrtReal(sqr_min_max.first / sqr_min_max.second);
}


real_t Polyhedron::computeTetrSimpleSqrQuality(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return sqr_min_max.first / sqr_min_max.second;
}




front::Edge* Polyhedron::currentFrontEdge(real_t maxCompl) const
{
    real_t cur_max_compl = static_cast<real_t>(0.0);
    front::Edge* cur_max_f_edge = nullptr;
    for (auto& f_edge : m_frontEdges)
    {
        real_t cur_compl = f_edge->complexity();
        if (cur_compl > cur_max_compl &&
            cur_compl < maxCompl)
        {
            cur_max_compl = cur_compl;
            cur_max_f_edge = f_edge;
        }
    }

    return cur_max_f_edge;
}


bool Polyhedron::exhaustWithoutNewVertPriorityPredicate(front::Edge* curFEdge)
{
    if (curFEdge->angleExCos() > cosDeg<70, real_t>)
        return true;

    auto opp_verts = curFEdge->findOppVerts();
    if (curFEdge->findOppEdge() ||
        (curFEdge->angleExCos() < cosDeg<70, real_t> &&
         curFEdge->angleExCos() > cosDeg<100, real_t> &&
         (*std::get<1>(opp_verts) - *std::get<0>(opp_verts)).sqrMagnitude() <= m_preferredLength * m_preferredLength))
        return true;

    if (frontSplitCheck(curFEdge))
        return true;

    return false;
}


bool Polyhedron::exhaustWithNewVertPriorityPredicate(front::Edge* currentFrontEdge)
{
    if (currentFrontEdge->angleExCos() < cosDeg<120, real_t>)
        return true;

    return false;
}

 
Polyhedron::ExhaustType Polyhedron::computeExhaustionTypeQualityPriority(
    front::Edge* currentFrontEdge,
    front::Face*& out_withNWFrontFace, Vec*& out_withNWNewVertPos)
{
    if (frontSplitCheck(currentFrontEdge))
        return ExhaustType::WithoutNewVert;
    
    if (parallelFacesCheck(currentFrontEdge) ||
        edgeIntersectionCheck(currentFrontEdge) ||
        facesIntersectionCheck(currentFrontEdge) ||
        anyVertInsidePotentialTetrCheck(currentFrontEdge))
    {
        return ExhaustType::WithNewVert;
    }

    auto opp_verts = currentFrontEdge->findOppVerts();
    real_t without_nv_quality = computeTetrSimpleSqrQuality(
        currentFrontEdge->edge->verts[0]->pos(),
        currentFrontEdge->edge->verts[1]->pos(),
        std::get<0>(opp_verts)->pos(),
        std::get<1>(opp_verts)->pos());

    front::Face* fface = chooseFaceForExhaustionWithNewVert(currentFrontEdge);
    Vec new_vert_pos;
    if (!tryComputeNewVertPos(fface, new_vert_pos))
        return ExhaustType::DontExhaust;

    real_t with_nv_quality = computeTetrSimpleSqrQuality(
        fface->face->edges[0]->verts[0]->pos(),
        fface->face->edges[0]->verts[1]->pos(),
        fface->face->findVertNot(fface->face->edges[0])->pos(),
        new_vert_pos);

    if (without_nv_quality > with_nv_quality)
        return ExhaustType::WithoutNewVert;

    out_withNWFrontFace = fface;
    out_withNWNewVertPos = new Vec(new_vert_pos);
    return ExhaustType::WithNewVert;
}




Vec Polyhedron::computeNormalInTetr(const front::Face* fFace, const Vec& oppVertPos) const
{
    return computeNormalInTetr(fFace->face->edges[0]->verts[0]->pos(),
                               fFace->face->edges[0]->verts[1]->pos(),
                               fFace->face->findVertNot(fFace->face->edges[0])->pos(),
                               oppVertPos);
}


Vec Polyhedron::computeNormalInTetr(const front::Face* fFace, const pmg::Edge* oneOfRemainingEdges) const
{
    Vec opp_pos = fFace->face->contains(oneOfRemainingEdges->verts[0]) ?
        oneOfRemainingEdges->verts[1]->pos() :
        oneOfRemainingEdges->verts[0]->pos();

    return computeNormalInTetr(fFace, opp_pos);
}


Vec Polyhedron::computeNormalInTetr(const Vec& fFacePos0, const Vec& fFacePos1, const Vec& fFacePos2, const Vec& oppVertPos) const
{
    Vec normal = Vec::cross(
        fFacePos0 - fFacePos2,
        fFacePos1 - fFacePos2).normalize();
    if (Vec::dot(normal, oppVertPos - fFacePos2) > static_cast<real_t>(0.0))
        normal *= static_cast<real_t>(-1.0);

    return normal;
}




void Polyhedron::setFEdgesInFrontSplit(const front::Edge* fEdge, front::Edge* newOppFEdges[2], front::Face* newFFaces[2], pair_ff oppFFaces) const
{
//    real_t cpa_times[2][2];
//    cpa_times[0][0] = tva::spatalgs::cpaTime(newFFaces[0]->computeCenter(), newFFaces[0]->normal,
//                                             oppFFaces.first->computeCenter(), oppFFaces.first->normal);
//    cpa_times[0][1] = tva::spatalgs::cpaTime(newFFaces[0]->computeCenter(), newFFaces[0]->normal,
//                                             oppFFaces.second->computeCenter(), oppFFaces.second->normal);
//    cpa_times[1][0] = tva::spatalgs::cpaTime(newFFaces[1]->computeCenter(), newFFaces[1]->normal,
//                                             oppFFaces.first->computeCenter(), oppFFaces.first->normal);
//    cpa_times[1][1] = tva::spatalgs::cpaTime(newFFaces[1]->computeCenter(), newFFaces[1]->normal,
//                                             oppFFaces.second->computeCenter(), oppFFaces.second->normal);

//    std::pair<front::Edge*, front::Face*> pair0;
//    std::pair<front::Edge*, front::Face*> pair1;
//    bool trivial_solution = false;

//    if (cpa_times[0][0] > 1e-6 && cpa_times[1][1] > 1e-6)
//    {
//        pair0 = { newOppFEdges[0], oppFFaces.first };
//        pair1 = { newOppFEdges[1], oppFFaces.second };
//        trivial_solution = true;
//    }

//    if (cpa_times[0][1] > 1e-6 && cpa_times[1][0] > 1e-6)
//    {
//        pair0 = { newOppFEdges[0], oppFFaces.second };
//        pair1 = { newOppFEdges[1], oppFFaces.first };
//        trivial_solution = true;
//    }

//    if (trivial_solution)
//    {
//        pair0.first->fillAdjFFaces(newFFaces[0], pair0.second);
//        pair0.second->addFEdge(pair0.first);

//        pair1.first->fillAdjFFaces(newFFaces[1], pair1.second);
//        pair1.second->addFEdge(pair1.first);
//        return;
//    }

//    std::cout << "\nPolyhedron::setFEdgesInFrontSplit here is not trivial solution...";
//    std::cin.get();

    pmg::Edge* opp_edge = newOppFEdges[0]->edge;
    Vec opp_verts_poses[2];
    opp_verts_poses[0] = opp_edge->verts[0]->pos();
    opp_verts_poses[1] = opp_edge->verts[1]->pos();
    Vec main_vert0_pos  = newFFaces[0]->face->contains(fEdge->edge->verts[0]) ?
                fEdge->edge->verts[0]->pos() : fEdge->edge->verts[1]->pos();
    Vec main_vert1_pos  = newFFaces[1]->face->contains(fEdge->edge->verts[0]) ?
                fEdge->edge->verts[0]->pos() : fEdge->edge->verts[1]->pos();
    Vec main_vert0_proj = spatalgs::project(main_vert0_pos, opp_verts_poses[0], opp_verts_poses[1]);
    Vec main_vert1_proj = spatalgs::project(main_vert1_pos, opp_verts_poses[0], opp_verts_poses[1]);
    Vec main_vec0 = main_vert0_pos - main_vert0_proj;
    Vec main_vec1 = main_vert1_pos - main_vert1_proj;

    Vec adj_opp_pos0 = oppFFaces.first->face->findVertNot(opp_edge)->pos();
    Vec adj_opp_pos1 = oppFFaces.second->face->findVertNot(opp_edge)->pos();
    Vec adj_opp_proj0 = spatalgs::project(adj_opp_pos0, opp_verts_poses[0], opp_verts_poses[1]);
    Vec adj_opp_proj1 = spatalgs::project(adj_opp_pos1, opp_verts_poses[0], opp_verts_poses[1]);
    Vec adj_vec0 = adj_opp_pos0 - adj_opp_proj0;
    Vec adj_vec1 = adj_opp_pos1 - adj_opp_proj1;

    real_t coses[2][2];
    coses[0][0] = Vec::cos(main_vec0, adj_vec0);
    coses[0][1] = Vec::cos(main_vec0, adj_vec1);
    coses[1][0] = Vec::cos(main_vec1, adj_vec0);
    coses[1][1] = Vec::cos(main_vec1, adj_vec1);

    if (coses[0][0] > coses[0][1] && coses[1][1] > coses[1][0])
    {
        newOppFEdges[0]->fillAdjFFaces(newFFaces[0], oppFFaces.first);
        oppFFaces.first->addFEdge(newOppFEdges[0]);

        newOppFEdges[1]->fillAdjFFaces(newFFaces[1], oppFFaces.second);
        oppFFaces.second->addFEdge(newOppFEdges[1]);
        return;
    }

    if (coses[0][1] >= coses[0][0] && coses[1][0] >= coses[1][1])
    {
        newOppFEdges[0]->fillAdjFFaces(newFFaces[0], oppFFaces.second);
        oppFFaces.second->addFEdge(newOppFEdges[0]);

        newOppFEdges[1]->fillAdjFFaces(newFFaces[1], oppFFaces.first);
        oppFFaces.first->addFEdge(newOppFEdges[1]);
        return;
    }

    int best_ofi = -1;
    if (coses[0][0] > coses[0][1] && coses[1][0] > coses[1][1])
        best_ofi = 0;

    if (coses[0][1] >= coses[0][0] && coses[1][1] >= coses[1][0])
        best_ofi = 1;

    int worst_ofi = best_ofi == 0 ? 1 : 0;

    front::Face* opp_ffaces[2] = { oppFFaces.first, oppFFaces.second };
    if (coses[0][best_ofi] > coses[1][best_ofi])
    {
        newOppFEdges[0]->fillAdjFFaces(newFFaces[0], opp_ffaces[best_ofi]);
        opp_ffaces[best_ofi]->addFEdge(newOppFEdges[0]);

        newOppFEdges[1]->fillAdjFFaces(newFFaces[1], opp_ffaces[worst_ofi]);
        opp_ffaces[worst_ofi]->addFEdge(newOppFEdges[1]);
    }
    else
    {
        newOppFEdges[0]->fillAdjFFaces(newFFaces[0], opp_ffaces[worst_ofi]);
        opp_ffaces[worst_ofi]->addFEdge(newOppFEdges[0]);

        newOppFEdges[1]->fillAdjFFaces(newFFaces[1], opp_ffaces[best_ofi]);
        opp_ffaces[best_ofi]->addFEdge(newOppFEdges[1]);
    }

    return;
}


void Polyhedron::exhaustFrontCollapse(front::Edge *fEdge, front::Edge *oppFEdge)
{
    auto fedge_adj_faces = fEdge->getAdjFFaces();
    auto opp_ffaces   = oppFEdge->getAdjFFaces();

    m_innerTetrs.push_back(new Tetr(
        fEdge->edge->verts[0],
        fEdge->edge->verts[1],
        oppFEdge->edge->verts[0],
        oppFEdge->edge->verts[1]));

    auto new_tetr = new Tetr(fEdge->edge->verts[0],
                             fEdge->edge->verts[1],
                             oppFEdge->edge->verts[0],
                             oppFEdge->edge->verts[1]);
    m_innerTetrs.push_back(new_tetr);

//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(PolyhedralSet::FileType::Obj, "debug.obj");


    std::vector<front::Edge*> fedges_to_erase;
    fedges_to_erase.reserve(4);
    for (auto& fedge : opp_ffaces.first->fEdges)
        if (fedge->edge->contains(fEdge->edge->verts[0]) || fedge->edge->contains(fEdge->edge->verts[1]))
            fedges_to_erase.push_back(fedge);

    for (auto& fedge : opp_ffaces.second->fEdges)
        if (fedge->edge->contains(fEdge->edge->verts[0]) || fedge->edge->contains(fEdge->edge->verts[1]))
            fedges_to_erase.push_back(fedge);


    removeFromFront(fedge_adj_faces.first);
    removeFromFront(fedge_adj_faces.second);
    removeFromFront(opp_ffaces.first);
    removeFromFront(opp_ffaces.second);
    removeFromFront(fEdge);
    removeFromFront(oppFEdge);
    for (auto& fedge : fedges_to_erase)
        removeFromFront(fedge);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustFrontCollapse";
#endif
}


void Polyhedron::exhaustFrontSplit(front::Edge* fEdge, front::Edge* oppFEdge)
{
    // This volatile helps to avoid computational error.
    volatile auto adj_ffaces = fEdge->getAdjFFaces();

    Vert* opp_verts[2];
    opp_verts[0] = adj_ffaces.first->face->findVertNot(fEdge->edge);
    opp_verts[1] = adj_ffaces.second->face->findVertNot(fEdge->edge);

    auto opp_edge    = oppFEdge->edge;
    auto opp_ffaces = oppFEdge->getAdjFFaces();

    opp_ffaces.first->removeFEdge(oppFEdge);
    opp_ffaces.second->removeFEdge(oppFEdge);
    removeFromFront(oppFEdge);
    front::Edge* new_opp_fedges[2];
    new_opp_fedges[0] = addToFront(opp_edge, false);
    new_opp_fedges[1] = addToFront(opp_edge, false);

    front::Edge* new_tetr_fedges[3];
    new_tetr_fedges[0] = nullptr;
    new_tetr_fedges[1] = nullptr;
    new_tetr_fedges[2] = new_opp_fedges[0];

    front::Face* new_ffaces[2];

    for (auto& fedge : adj_ffaces.first->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[0]) &&
            (fedge->edge->contains(opp_edge->verts[0]) ||
             fedge->edge->contains(opp_edge->verts[1])))
        {
            new_tetr_fedges[0] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFace(adj_ffaces.first);
            fedge->removeAdjFFace(adj_ffaces.second);
            break;
        }
    }
    for (auto& fedge : adj_ffaces.second->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[0]) &&
            (fedge->edge->contains(opp_edge->verts[0]) ||
             fedge->edge->contains(opp_edge->verts[1])))
        {
            new_tetr_fedges[1] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFace(adj_ffaces.first);
            fedge->removeAdjFFace(adj_ffaces.second);
            break;
        }
    }
    new_ffaces[0] = addToFront(new Face(new_tetr_fedges[0]->edge,
                                          new_tetr_fedges[1]->edge,
                                          new_tetr_fedges[2]->edge));
    new_ffaces[0]->addFEdge(new_tetr_fedges[0]);
    new_ffaces[0]->addFEdge(new_tetr_fedges[1]);
    new_ffaces[0]->addFEdge(new_tetr_fedges[2]);
    new_ffaces[0]->normal = computeNormalInTetr(new_ffaces[0], fEdge->edge);


    new_tetr_fedges[2] = new_opp_fedges[1];

    for (auto& fedge : adj_ffaces.first->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[1]) &&
            (fedge->edge->contains(opp_edge->verts[0]) ||
             fedge->edge->contains(opp_edge->verts[1])))
        {
            new_tetr_fedges[0] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFace(adj_ffaces.first);
            fedge->removeAdjFFace(adj_ffaces.second);
            break;
        }
    }
    for (auto& fedge : adj_ffaces.second->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[1]) &&
            (fedge->edge->contains(opp_edge->verts[0]) ||
             fedge->edge->contains(opp_edge->verts[1])))
        {
            new_tetr_fedges[1] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFace(adj_ffaces.first);
            fedge->removeAdjFFace(adj_ffaces.second);
            break;
        }
    }
    new_ffaces[1] = addToFront(new Face(new_tetr_fedges[0]->edge,
                                          new_tetr_fedges[1]->edge,
                                          new_tetr_fedges[2]->edge));
    new_ffaces[1]->addFEdge(new_tetr_fedges[0]);
    new_ffaces[1]->addFEdge(new_tetr_fedges[1]);
    new_ffaces[1]->addFEdge(new_tetr_fedges[2]);
    new_ffaces[1]->normal = computeNormalInTetr(new_ffaces[1], fEdge->edge);


    setFEdgesInFrontSplit(fEdge, new_opp_fedges, new_ffaces, opp_ffaces);

//    m_innerTetrs.push_back(new Tetr(
//        fEdge->edge->verts[0],
//        fEdge->edge->verts[1],
//        opp_verts[0],
//        opp_verts[1]));

    auto new_tetr = new Tetr(fEdge->edge->verts[0],
                             fEdge->edge->verts[1],
                             opp_verts[0],
                             opp_verts[1]);
    m_innerTetrs.push_back(new_tetr);

//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(PolyhedralSet::FileType::Obj, "debug.obj");


    removeFromFront(adj_ffaces.first);
    removeFromFront(adj_ffaces.second);
    removeFromFront(fEdge);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustFrontSplit";
#endif
}


void Polyhedron::exhaustWithoutNewVertOppEdgeExists(front::Edge* fEdge, front::Edge* oppFEdge)
{
    auto opp_ffaces = oppFEdge->getAdjFFaces();

    front::Face* main_ffaces[3];
    Vert* main_vert;
    auto fedge_adj_faces = fEdge->getAdjFFaces();
    main_ffaces[0] = std::get<0>(fedge_adj_faces);
    main_ffaces[1] = std::get<1>(fedge_adj_faces);

    if (std::get<0>(opp_ffaces)->face->contains(fEdge->edge->verts[0]))
    {
        main_ffaces[2] = std::get<0>(opp_ffaces);
        main_vert = fEdge->edge->verts[0];
    }
    else if (std::get<0>(opp_ffaces)->face->contains(fEdge->edge->verts[1]))
    {
        main_ffaces[2] = std::get<0>(opp_ffaces);
        main_vert = fEdge->edge->verts[1];
    }
    else if (std::get<1>(opp_ffaces)->face->contains(fEdge->edge->verts[0]))
    {
        main_ffaces[2] = std::get<1>(opp_ffaces);
        main_vert = fEdge->edge->verts[0];
    }
    else
    {
        main_ffaces[2] = std::get<1>(opp_ffaces);
        main_vert = fEdge->edge->verts[1];
    }

    front::Edge* new_tetr_fedges[3];
    for (int i = 0; i < 3; i++)
    {
        new_tetr_fedges[i] = main_ffaces[i]->findFEdgeNot(main_vert);
        new_tetr_fedges[i]->refreshAngleData();
        new_tetr_fedges[i]->removeAdjFFace(main_ffaces[i]);
    }

    auto new_fface = addToFront(new Face(new_tetr_fedges[0]->edge,
                                           new_tetr_fedges[1]->edge,
                                           new_tetr_fedges[2]->edge));
    new_fface->addFEdge(new_tetr_fedges[0]);
    new_fface->addFEdge(new_tetr_fedges[1]);
    new_fface->addFEdge(new_tetr_fedges[2]);
    new_fface->normal = computeNormalInTetr(new_fface, main_vert->pos());

//    m_innerTetrs.push_back(new Tetr(
//        fEdge->edge->verts[0],
//        fEdge->edge->verts[1],
//        oppFEdge->edge->verts[0],
//        oppFEdge->edge->verts[1]));

    auto new_tetr = new Tetr(fEdge->edge->verts[0],
                             fEdge->edge->verts[1],
                             oppFEdge->edge->verts[0],
                             oppFEdge->edge->verts[1]);
    m_innerTetrs.push_back(new_tetr);

//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(PolyhedralSet::FileType::Obj, "debug.obj");

    std::vector<front::Edge*> erased_fedges;
    erased_fedges.reserve(3);
    for (auto& fface : main_ffaces)
    {
        for (auto& fedge : fface->fEdges)
        {
            if (std::find(erased_fedges.begin(), erased_fedges.end(), fedge) == erased_fedges.end() &&
                fedge->edge->contains(main_vert))
            {
                removeFromFront(fedge);
                erased_fedges.push_back(fedge);
            }
        }
    }

    for (auto& fface : main_ffaces)
        removeFromFront(fface);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustWithoutNewVertOppEdgeExists";
#endif
}


void Polyhedron::exhaustWithoutNewVertOppEdgeDontExists(front::Edge* fEdge)
{
    // This volatile helps to avoid computational error.
    volatile auto adj_ffaces = fEdge->getAdjFFaces();

    Vert* opp_verts[2];
    opp_verts[0] = adj_ffaces.first->face->findVertNot(fEdge->edge);
    opp_verts[1] = adj_ffaces.second->face->findVertNot(fEdge->edge);
    
    front::Edge* opp_fedge = addToFront(new Edge(opp_verts[0], opp_verts[1]));

    front::Edge* new_tetr_fedges[3];
    new_tetr_fedges[0] = nullptr;
    new_tetr_fedges[1] = nullptr;
    new_tetr_fedges[2] = opp_fedge;

    for (auto& fedge : adj_ffaces.first->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[0]) &&
            (fedge->edge->contains(opp_fedge->edge->verts[0]) ||
             fedge->edge->contains(opp_fedge->edge->verts[1])))
        {
            new_tetr_fedges[0] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFace(adj_ffaces.first);
//            fedge->removeAdjFFace(adj_ffaces.second);
            break;
        }
    }
    for (auto& fedge : adj_ffaces.second->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[0]) &&
            (fedge->edge->contains(opp_fedge->edge->verts[0]) ||
             fedge->edge->contains(opp_fedge->edge->verts[1])))
        {
            new_tetr_fedges[1] = fedge;
            fedge->refreshAngleData();
//            fedge->removeAdjFFace(adj_ffaces.first);
            fedge->removeAdjFFace(adj_ffaces.second);
            break;
        }
    }
    auto new_fface = addToFront(new Face(new_tetr_fedges[0]->edge,
                                           new_tetr_fedges[1]->edge,
                                           new_tetr_fedges[2]->edge));
    new_fface->addFEdge(new_tetr_fedges[0]);
    new_fface->addFEdge(new_tetr_fedges[1]);
    new_fface->addFEdge(new_tetr_fedges[2]);
    new_fface->normal = computeNormalInTetr(new_fface, fEdge->edge);

    for (auto& fedge : adj_ffaces.first->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[1]) &&
            (fedge->edge->contains(opp_fedge->edge->verts[0]) ||
             fedge->edge->contains(opp_fedge->edge->verts[1])))
        {
            new_tetr_fedges[0] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFace(adj_ffaces.first);
            fedge->removeAdjFFace(adj_ffaces.second);
            break;
        }
    }
    for (auto& fedge : adj_ffaces.second->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[1]) &&
            (fedge->edge->contains(opp_fedge->edge->verts[0]) ||
             fedge->edge->contains(opp_fedge->edge->verts[1])))
        {
            new_tetr_fedges[1] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFace(adj_ffaces.first);
            fedge->removeAdjFFace(adj_ffaces.second);
            break;
        }
    }
    new_fface = addToFront(new Face(new_tetr_fedges[0]->edge,
                                      new_tetr_fedges[1]->edge,
                                      new_tetr_fedges[2]->edge));
    new_fface->addFEdge(new_tetr_fedges[0]);
    new_fface->addFEdge(new_tetr_fedges[1]);
    new_fface->addFEdge(new_tetr_fedges[2]);
    new_fface->normal = computeNormalInTetr(new_fface, fEdge->edge);

//    m_innerTetrs.push_back(new Tetr(
//        fEdge->edge->verts[0],
//        fEdge->edge->verts[1],
//        opp_verts[0],
//        opp_verts[1]));

    auto new_tetr = new Tetr(fEdge->edge->verts[0],
                             fEdge->edge->verts[1],
                             opp_verts[0],
                             opp_verts[1]);
    m_innerTetrs.push_back(new_tetr);

//    std::cout << std::endl << fEdge->edge->magnitude() << ' ' << fEdge->computeAngle() * 180.0 / M_PI;
//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(PolyhedralSet::FileType::Obj, "debug.obj");

    
    removeFromFront(adj_ffaces.first);
    removeFromFront(adj_ffaces.second);
    removeFromFront(fEdge);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustWithoutNewVertOppEdgeDontExists";
#endif
}


void Polyhedron::exhaustWithoutNewVert(front::Edge* fEdge, bool oppEdgeExistence, front::Edge* oppFEdge)
{
    front::Edge* opp_fedge = nullptr;
    if (oppEdgeExistence && oppFEdge)
        opp_fedge = oppFEdge;
    else if (oppEdgeExistence)
        opp_fedge = fEdge->findOppEdge();

    if (/*(oppFEdge && oppEdgeExistence) ||*/
        opp_fedge)
    {
        if (frontCollapseCheck(fEdge, opp_fedge))
        {
            exhaustFrontCollapse(fEdge, opp_fedge);
        }
        else if (frontSplitCheck(fEdge, opp_fedge))
        {
            exhaustFrontSplit(fEdge, opp_fedge);
        }
        else
        {
            exhaustWithoutNewVertOppEdgeExists(fEdge, opp_fedge);
        }
    }
    else
    {
        exhaustWithoutNewVertOppEdgeDontExists(fEdge);
    }
}




bool Polyhedron::tryComputeNewVertPosType3(front::Face* fFace, Vec& out_pos)
{
    front::Edge* main_fedges[2];
    main_fedges[0] = fFace->fEdges[0];
    main_fedges[1] = fFace->fEdges[1];
    Edge* main_edges[2];
    main_edges[0] = main_fedges[0]->edge;
    main_edges[1] = main_fedges[1]->edge;

    auto main_vert = main_edges[0]->verts[0] == main_edges[1]->verts[0] ?
        main_edges[0]->verts[0] :
        (main_edges[0]->verts[0] == main_edges[1]->verts[1] ?
            main_edges[0]->verts[0] :
            main_edges[0]->verts[1]);
    auto third_edge = fFace->face->findEdgeNot(main_vert);
    auto third_f_edge = fFace->findFEdge(third_edge);

    auto v0 = main_edges[0]->verts[0] == main_vert ? main_edges[0]->verts[1] : main_edges[0]->verts[0];
    auto v1 = main_edges[1]->verts[0] == main_vert ? main_edges[1]->verts[1] : main_edges[1]->verts[0];
    auto v2 = main_vert;

    auto adj_ffaces = main_fedges[0]->getAdjFFaces();
    auto fn0 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    adj_ffaces = main_fedges[1]->getAdjFFaces();
    auto fn1 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    adj_ffaces = third_f_edge->getAdjFFaces();
    auto fn2 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);

    const Vec& v0_pos = v0->pos();
    const Vec& v1_pos = v1->pos();
    const Vec& v2_pos = v2->pos();

    Vec n_m  = fFace->normal;
    Vec n_n0 = fn0->normal;
    Vec n_n1 = fn1->normal;
    Vec n_n2 = fn2->normal;

    Vec e_mn0 = (n_m + n_n0).normalize();
    Vec e_mn1 = (n_m + n_n1).normalize();
    Vec e_mn2 = n_m + n_n2;

    Vec np_mn0 = Vec::cross(v2_pos - v0_pos, e_mn0);
    Vec np_mn1 = Vec::cross(v2_pos - v1_pos, e_mn1);

    Vec e = Vec::cross(np_mn0, np_mn1);

    Vec new_pos = spatalgs::lineIntersectPlane(v2_pos, e, v0_pos, v1_pos, v0_pos + e_mn2);

    Vec v0_to_np = new_pos - v0_pos;
    Vec v1_to_np = new_pos - v1_pos;
    Vec v2_to_np = new_pos - v2_pos;
    real_t sum_magns[4];
    int i = 0;
    for (auto& fface : { fn0, fn1, fn2, fFace })
        sum_magns[i++] =  fface->face->edges[0]->magnitude()
                        + fface->face->edges[1]->magnitude()
                        + fface->face->edges[2]->magnitude();
    real_t av_magn = (sum_magns[0] + sum_magns[1] + sum_magns[2] + sum_magns[3]) / static_cast<real_t>(12.0);
    if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos /*+ NOT_TOO_CLOSE * v0_to_np*/) ||
        segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos /*+ NOT_TOO_CLOSE * v1_to_np*/) ||
        segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos /*+ NOT_TOO_CLOSE * v2_to_np*/) ||
        doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * av_magn) ||
//        !vertInsideFrontCheck(new_pos) ||
        faceIntersectionCheck(v0, v1, new_pos) ||
        faceIntersectionCheck(v0, v2, new_pos) ||
        faceIntersectionCheck(v1, v2, new_pos))
        return false;

//    std::cout << std::endl << "Type3";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::tryComputeNewVertPosType2(front::Face* fFace, Vec& out_pos, int smallAngleIndex0, int smallAngleIndex1)
{
    front::Edge* main_fedges[2];
    main_fedges[0] = fFace->fEdges[smallAngleIndex0];
    main_fedges[1] = fFace->fEdges[smallAngleIndex1];
    Edge* main_edges[2];
    main_edges[0] = main_fedges[0]->edge;
    main_edges[1] = main_fedges[1]->edge;
    auto main_vert = main_edges[0]->verts[0] == main_edges[1]->verts[0] ?
        main_edges[0]->verts[0] :
        (main_edges[0]->verts[0] == main_edges[1]->verts[1] ?
            main_edges[0]->verts[0] :
            main_edges[0]->verts[1]);

    auto v0 = main_edges[0]->verts[0] == main_vert ? main_edges[0]->verts[1] : main_edges[0]->verts[0];
    auto v1 = main_edges[1]->verts[0] == main_vert ? main_edges[1]->verts[1] : main_edges[1]->verts[0];
    auto v2 = main_vert;

    auto adj_ffaces = main_fedges[0]->getAdjFFaces();
    auto fn0 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn0 = fn0->face->findVertNot(main_fedges[0]->edge);
    adj_ffaces = main_fedges[1]->getAdjFFaces();
    auto fn1 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn1 = fn1->face->findVertNot(main_fedges[1]->edge);

    Vec v0_pos = v0->pos();
    Vec v1_pos = v1->pos();
    Vec v2_pos = v2->pos();
    Vec vn0_pos = vn0->pos();
    Vec vn1_pos = vn1->pos();

    Vec n_m = fFace->normal;
    Vec n_n0 = fn0->normal;
    Vec n_n1 = fn1->normal;
    
    Vec e_mn0 = (n_m + n_n0).normalize();
    Vec e_mn1 = (n_m + n_n1).normalize();

    Vec np_mn0 = Vec::cross(v2_pos - v0_pos, e_mn0);
    Vec np_mn1 = Vec::cross(v2_pos - v1_pos, e_mn1);

    Vec e = Vec::cross(np_mn0, np_mn1).normalize();
    if (Vec::dot(e, n_m) < static_cast<real_t>(0.0)) e *= static_cast<real_t>(-1.0);

//    real_t fface_area = fFace->Face->computeArea();
//    real_t sp = SQRT3_2 * 0.5 * m_preferredLength * m_preferredLength;
//    real_t sf0 = 0.5 * (fface_area + fn0->Face->computeArea());
//    real_t sf1 = 0.5 * (fface_area + fn1->Face->computeArea());
//    real_t raw_deform0 = K_D * (sp - sf0);
//    real_t raw_deform1 = K_D * (sp - sf1);
//    real_t deform0 = raw_deform0 < sf0 * K_MAXD ? raw_deform0 : sf0 * K_MAXD;
//    real_t deform1 = raw_deform1 < sf1 * K_MAXD ? raw_deform1 : sf1 * K_MAXD;
//    real_t sc0 = sf0 + deform0;
//    real_t sc1 = sf1 + deform1;
//    real_t x0_2 = sc0 / Vec::cross(v2->pos() - v0->pos(), e).magnitude();
//    real_t x1_2 = sc1 / Vec::cross(v2->pos() - v1->pos(), e).magnitude();
//    Vec new_pos = v2_pos + (x0_2 + x1_2) * e;
    real_t lm0 = main_edges[0]->magnitude();
    real_t lm1 = main_edges[1]->magnitude();
    real_t l0 = (v1_pos - v0_pos).magnitude();
    real_t l1 = (vn0_pos - v0_pos).magnitude();
    real_t l2 = (vn0_pos - v2_pos).magnitude();
    real_t l3 = (vn1_pos - v1_pos).magnitude();
    real_t l4 = (vn1_pos - v2_pos).magnitude();
    real_t av_magn = (lm0 + lm1 + l0 + l1 + l2 + l3 + l4) / static_cast<real_t>(7.0);
    real_t raw_deform = K_D * (m_preferredLength - av_magn);
    real_t deform = raw_deform < av_magn * K_MAXD ? raw_deform : av_magn * K_MAXD;
    real_t magn_d = av_magn + deform;
    Vec new_pos = v2_pos + magn_d * e;

    Vec v0_to_np = new_pos - v0_pos;
    Vec v1_to_np = new_pos - v1_pos;
    Vec v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos /*+ NOT_TOO_CLOSE * v0_to_np*/) ||
        segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos /*+ NOT_TOO_CLOSE * v1_to_np*/) ||
        segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos /*+ NOT_TOO_CLOSE * v2_to_np*/) ||
        doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * magn_d) ||
        !vertInsideFrontCheck(new_pos) ||
        faceIntersectionCheck(v0, v1, new_pos) ||
        faceIntersectionCheck(v0, v2, new_pos) ||
        faceIntersectionCheck(v1, v2, new_pos))
    {
//        x0_2 = sf0 / Vec::cross(v2->pos() - v0->pos(), e).magnitude();
//        x1_2 = sf1 / Vec::cross(v2->pos() - v1->pos(), e).magnitude();
        new_pos = new_pos = v2_pos + av_magn * e;
        v0_to_np = new_pos - v0_pos;
        v1_to_np = new_pos - v1_pos;
        v2_to_np = new_pos - v2_pos;
        if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos /*+ NOT_TOO_CLOSE * v0_to_np*/) ||
            segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos /*+ NOT_TOO_CLOSE * v1_to_np*/) ||
            segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos /*+ NOT_TOO_CLOSE * v2_to_np*/) ||
            doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * av_magn) ||
//            !vertInsideFrontCheck(new_pos) ||
            faceIntersectionCheck(v0, v1, new_pos) ||
            faceIntersectionCheck(v0, v2, new_pos) ||
            faceIntersectionCheck(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type2";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::tryComputeNewVertPosType1(front::Face* fFace, Vec& out_pos, int smallAngleIndex)
{
    auto main_f_edge = fFace->fEdges[smallAngleIndex];
    auto main_edge   = fFace->fEdges[smallAngleIndex]->edge;
    Vec main_edge_poses[2];
    main_edge_poses[0] = main_edge->verts[0]->pos();
    main_edge_poses[1] = main_edge->verts[1]->pos();

    auto v0 = main_edge->verts[0];
    auto v1 = main_edge->verts[1];
    auto v2 = fFace->face->findVertNot(main_edge);

    auto adj_ffaces = main_f_edge->getAdjFFaces();
    auto fn = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn = fn->face->findVertNot(main_edge);

    Vec v0_pos = v0->pos();
    Vec v1_pos = v1->pos();
    Vec v2_pos = v2->pos();
    Vec vn_pos = vn->pos();

    Vec v2pr = spatalgs::project(v2_pos, v0_pos, v1_pos);
    Vec vnpr = spatalgs::project(vn_pos, v0_pos, v1_pos);

    Vec c = static_cast<real_t>(0.25) * (v0_pos + v1_pos + v2pr + vnpr);
    Vec e = (fFace->normal + fn->normal).normalize();
    real_t me_magn = main_edge->magnitude();
    real_t l0 = (v2_pos - v0_pos).magnitude();
    real_t l1 = (v2_pos - v1_pos).magnitude();
    real_t l2 = (vn_pos - v0_pos).magnitude();
    real_t l3 = (vn_pos - v1_pos).magnitude();
    real_t av_magn = static_cast<real_t>(0.2) * (me_magn + l0 + l1 + l2 + l3);
    real_t v0_c_dist = (c - v0_pos).magnitude();
    real_t raw_deform = K_D * (m_preferredLength - av_magn);
    real_t deform = raw_deform < av_magn * K_MAXD ? raw_deform : av_magn * K_MAXD;
    real_t magn_d = av_magn + deform;
    Vec new_pos = c + sqrtReal(magn_d * magn_d - v0_c_dist * v0_c_dist) * e;

    Vec v0_to_np = new_pos - v0_pos;
    Vec v1_to_np = new_pos - v1_pos;
    Vec v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos /*+ NOT_TOO_CLOSE * v0_to_np*/) ||
        segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos /*+ NOT_TOO_CLOSE * v1_to_np*/) ||
        segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos /*+ NOT_TOO_CLOSE * v2_to_np*/) ||
        doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * magn_d) ||
        !vertInsideFrontCheck(new_pos) ||
        faceIntersectionCheck(v0, v1, new_pos) ||
        faceIntersectionCheck(v0, v2, new_pos) ||
        faceIntersectionCheck(v1, v2, new_pos))
    {
        new_pos = c + sqrtReal(av_magn * av_magn - v0_c_dist * v0_c_dist) * e;
        v0_to_np = new_pos - v0_pos;
        v1_to_np = new_pos - v1_pos;
        v2_to_np = new_pos - v2_pos;
        if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos /*+ NOT_TOO_CLOSE * v0_to_np*/) ||
            segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos /*+ NOT_TOO_CLOSE * v1_to_np*/) ||
            segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos /*+ NOT_TOO_CLOSE * v2_to_np*/) ||
            doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * av_magn) ||
//            !vertInsideFrontCheck(new_pos) ||
            faceIntersectionCheck(v0, v1, new_pos) ||
            faceIntersectionCheck(v0, v2, new_pos) ||
            faceIntersectionCheck(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type1";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::tryComputeNewVertPosType0(front::Face* fFace, Vec& out_pos)
{
    real_t av_magn = ONE_3 * (
          fFace->face->edges[0]->magnitude()
        + fFace->face->edges[1]->magnitude()
        + fFace->face->edges[2]->magnitude());
    real_t raw_deform = K_D * (m_preferredLength - av_magn);
    real_t deform = raw_deform < av_magn * K_MAXD ? raw_deform : av_magn * K_MAXD;
    real_t magn_d = av_magn + deform;
    Vec new_pos = fFace->computeCenter() + sqrtReal(magn_d * magn_d - ONE_3 * av_magn * av_magn) * fFace->normal;

    auto v0 = fFace->face->edges[0]->verts[0];
    auto v1 = fFace->face->edges[0]->verts[1];
    auto v2 = fFace->face->findVertNot(fFace->face->edges[0]);
    const Vec& v0_pos = v0->pos();
    const Vec& v1_pos = v1->pos();
    const Vec& v2_pos = v2->pos();
    Vec v0_to_np = new_pos - v0_pos;
    Vec v1_to_np = new_pos - v1_pos;
    Vec v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos /*+ NOT_TOO_CLOSE * v0_to_np*/) ||
        segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos /*+ NOT_TOO_CLOSE * v1_to_np*/) ||
        segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos /*+ NOT_TOO_CLOSE * v2_to_np*/) ||
        doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * magn_d) ||
//        !vertInsideFrontCheck(new_pos) ||
        faceIntersectionCheck(v0, v1, new_pos) ||
        faceIntersectionCheck(v0, v2, new_pos) ||
        faceIntersectionCheck(v1, v2, new_pos))
    {
        new_pos = fFace->computeCenter() + sqrtReal(av_magn * av_magn - ONE_3 * av_magn * av_magn) * fFace->normal;
        v0_to_np = new_pos - v0_pos;
        v1_to_np = new_pos - v1_pos;
        v2_to_np = new_pos - v2_pos;
        if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos /*+ NOT_TOO_CLOSE * v0_to_np*/) ||
            segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos /*+ NOT_TOO_CLOSE * v1_to_np*/) ||
            segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos /*+ NOT_TOO_CLOSE * v2_to_np*/) ||
            doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * av_magn) ||
            !vertInsideFrontCheck(new_pos) ||
            faceIntersectionCheck(v0, v1, new_pos) ||
            faceIntersectionCheck(v0, v2, new_pos) ||
            faceIntersectionCheck(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type0";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::tryComputeNewVertPos(front::Face* fFace, Vec& out_pos)
{
    real_t angs_coses[3]
    { 
        fFace->fEdges[0]->angleExCos(),
        fFace->fEdges[1]->angleExCos(),
        fFace->fEdges[2]->angleExCos()
    };
    int indexes[3];
    int small_angs_num = 0;
    if (angs_coses[0] > cosDeg<140, real_t>) indexes[small_angs_num++] = 0;
    if (angs_coses[1] > cosDeg<140, real_t>) indexes[small_angs_num++] = 1;
    if (angs_coses[2] > cosDeg<140, real_t>) indexes[small_angs_num++] = 2;

    switch (small_angs_num)
    {
    case 0: return tryComputeNewVertPosType0(fFace, out_pos);
    case 1: return tryComputeNewVertPosType1(fFace, out_pos, indexes[0]);
    case 2: return tryComputeNewVertPosType2(fFace, out_pos, indexes[0], indexes[1]);
    case 3: return tryComputeNewVertPosType3(fFace, out_pos);
    }

    return true;
}




real_t Polyhedron::sqr4FaceArea(const front::Face* fFace) const
{
    Vec vec0 = *fFace->face->edges[0]->verts[1] - *fFace->face->edges[0]->verts[0];
    Vec vec1 = *fFace->face->edges[1]->verts[1] - *fFace->face->edges[1]->verts[0];
    return Vec::cross(vec0, vec1).sqrMagnitude();
}


front::Face* Polyhedron::chooseFaceForExhaustionWithNewVert(front::Edge* fEdge)
{
    auto adj_ffaces = fEdge->getAdjFFaces();

    return sqr4FaceArea(std::get<0>(adj_ffaces)) < sqr4FaceArea(std::get<1>(adj_ffaces)) ?
//    return std::get<0>(adj_ffaces)->computeQuality() < std::get<1>(adj_ffaces)->computeQuality() ?
        std::get<0>(adj_ffaces) :
        std::get<1>(adj_ffaces);
}


void Polyhedron::exhaustWithNewVert(front::Face* fFace, const Vec& vertPos)
{
    Vert* new_vert = new Vert(vertPos);
    m_innerVerts.push_back(new_vert);

    front::Edge* new_tetr_fedges[6];
    new_tetr_fedges[0] = fFace->fEdges[0];
    auto far_vert = fFace->face->findVertNot(new_tetr_fedges[0]->edge);
    new_tetr_fedges[1] = fFace->findFEdge(new_tetr_fedges[0]->edge->verts[1], far_vert);
    new_tetr_fedges[2] = fFace->findFEdge(new_tetr_fedges[0]->edge->verts[0], far_vert);
    new_tetr_fedges[3] = addToFront(new Edge(new_tetr_fedges[0]->edge->verts[0], new_vert));
    new_tetr_fedges[4] = addToFront(new Edge(new_tetr_fedges[0]->edge->verts[1], new_vert));
    new_tetr_fedges[5] = addToFront(new Edge(far_vert, new_vert));

    new_tetr_fedges[0]->refreshAngleData();
    new_tetr_fedges[1]->refreshAngleData();
    new_tetr_fedges[2]->refreshAngleData();

    new_tetr_fedges[0]->removeAdjFFace(fFace);
    new_tetr_fedges[1]->removeAdjFFace(fFace);
    new_tetr_fedges[2]->removeAdjFFace(fFace);
    
    auto new_fface = addToFront(new Face(new_tetr_fedges[0]->edge,
                                           new_tetr_fedges[3]->edge,
                                           new_tetr_fedges[4]->edge));
    new_fface->addFEdge(new_tetr_fedges[0]);
    new_fface->addFEdge(new_tetr_fedges[3]);
    new_fface->addFEdge(new_tetr_fedges[4]);
    new_fface->normal = computeNormalInTetr(new_fface, far_vert->pos());

    new_fface = addToFront(new Face(new_tetr_fedges[2]->edge,
                                      new_tetr_fedges[3]->edge,
                                      new_tetr_fedges[5]->edge));
    new_fface->addFEdge(new_tetr_fedges[2]);
    new_fface->addFEdge(new_tetr_fedges[3]);
    new_fface->addFEdge(new_tetr_fedges[5]);
    new_fface->normal = computeNormalInTetr(new_fface, new_tetr_fedges[0]->edge->verts[1]->pos());

    new_fface = addToFront(new Face(new_tetr_fedges[1]->edge,
                                       new_tetr_fedges[4]->edge,
                                       new_tetr_fedges[5]->edge));
    new_fface->addFEdge(new_tetr_fedges[1]);
    new_fface->addFEdge(new_tetr_fedges[4]);
    new_fface->addFEdge(new_tetr_fedges[5]);
    new_fface->normal = computeNormalInTetr(new_fface, new_tetr_fedges[0]->edge->verts[0]->pos());

//    m_innerTetrs.push_back(new Tetr(
//        new_tetr_fedges[0]->edge->verts[0],
//        new_tetr_fedges[0]->edge->verts[1],
//        new_tetr_fedges[5]->edge->verts[0],
//        new_tetr_fedges[5]->edge->verts[1]));

    auto new_tetr = new Tetr(new_tetr_fedges[0]->edge->verts[0],
                             new_tetr_fedges[0]->edge->verts[1],
                             new_tetr_fedges[5]->edge->verts[0],
                             new_tetr_fedges[5]->edge->verts[1]);
    m_innerTetrs.push_back(new_tetr);
//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(PolyhedralSet::FileType::Obj, "debug.obj");

    removeFromFront(fFace);
}




bool Polyhedron::tryExhaustWithoutNewVert(front::Edge* fEdge, bool oppEdgeExistence, front::Edge* oppEdge)
{
    if (parallelFacesCheck(fEdge) ||
        edgeIntersectionCheck(fEdge) ||
        facesIntersectionCheck(fEdge) ||
        anyVertInsidePotentialTetrCheck(fEdge))
        return false;

    exhaustWithoutNewVert(fEdge, oppEdgeExistence, oppEdge);
    return true;
}


bool Polyhedron::tryExhaustWithNewVert(front::Edge* frontEdge)
{
    if (parallelFacesCheck(frontEdge))
        return false;

    auto exhaust_fface = chooseFaceForExhaustionWithNewVert(frontEdge);
    Vec new_vert_pos;
    if (!tryComputeNewVertPos(exhaust_fface, new_vert_pos))
        return false;

    exhaustWithNewVert(exhaust_fface, new_vert_pos);
    return true;
}




bool Polyhedron::isFrontExhausted()
{
    if (m_frontFaces.size() == 0 &&
        m_frontEdges.size() == 0)
        return true;

    if ((m_frontFaces.size() == 0 && m_frontEdges.size() > 0) ||
        (m_frontFaces.size() > 0 && m_frontEdges.size() == 0))
        throw std::logic_error("Error in Polyhedron::isFrontExhausted. Front wasn't correctly exhausted.");

    return false;
}


void Polyhedron::processAngles()
{
#ifdef DEV_DEBUG
    int debug_i = 0;
#endif
    real_t max_compl = std::numeric_limits<real_t>::max();
    for (front::Edge* cur_fedge = currentFrontEdge(max_compl);; cur_fedge = currentFrontEdge(max_compl))
    {
        if (!cur_fedge)
            throw std::logic_error("pmg::Polyhedron::currentFrontEdge returned nullptr");
        
        if (exhaustWithoutNewVertPriorityPredicate(cur_fedge))
        {
            if (!tryExhaustWithoutNewVert(cur_fedge))
            {
                max_compl = cur_fedge->complexity();
                continue;
            }
        }
        else if (exhaustWithNewVertPriorityPredicate(cur_fedge))
        {
            if (!tryExhaustWithNewVert(cur_fedge))
            {
                max_compl = cur_fedge->complexity();
                continue;
            }
        }
        else
        {
            front::Face* exhaust_from_fface = nullptr;
            Vec* new_vert_pos = nullptr;
            switch (computeExhaustionTypeQualityPriority(cur_fedge, exhaust_from_fface, new_vert_pos))
            {
            case ExhaustType::WithoutNewVert:
                exhaustWithoutNewVert(cur_fedge);
                break;

            case ExhaustType::WithNewVert:
                if (new_vert_pos)
                {
                    exhaustWithNewVert(exhaust_from_fface, *new_vert_pos);
                    delete new_vert_pos;
                }
                else
                {
                    if (!tryExhaustWithNewVert(cur_fedge))
                    {
                        max_compl = cur_fedge->complexity();
                        continue;
                    }
                }
                break;

            case ExhaustType::DontExhaust:
                max_compl = cur_fedge->complexity();
#ifdef DEV_DEBUG
                std::cout << std::endl << "DontExhaust";
#endif
                continue;
            }
        }
        max_compl = std::numeric_limits<real_t>::max();

#ifdef DEV_DEBUG
//        if (debug_i++ >= 1200)
//        m_polycr->output(PolyhedralSet::FileType::Obj, "debug.obj");
//        if (debug_i++ >= 0)
//            std::cout << std::endl << debug_i - 1;

//        debug();
#endif

        if (isFrontExhausted())
            return;
    }
}


void Polyhedron::debug()
{
    Vec accum;
    for (auto& vert : m_innerVerts)
        accum += vert->pos();

    for (auto& svert : m_shellVerts)
        accum += svert->attachedVert->pos();

    for (auto& sedge : m_shellEdges)
        for (auto& vert : sedge->innerVerts())
            accum += vert->pos();

    for (auto& sface : m_shellFaces)
        for (auto& vert : sface->innerVerts())
            accum += vert->pos();

    std::cout << "\n{ " << accum.coors[0] + accum.coors[1] + accum.coors[2] << " }";
}




void Polyhedron::generateMesh(real_t preferredLen)
{
    m_preferredLength = preferredLen;
    initializeFront();
    computeFrontNormals();
    processAngles();
//    if (globalIntersectionCheck())
//        throw std::logic_error("Intersection error.\npmg::Polyhedron::globalIntersectionCheck returned true.");
    smoothMesh(20);
}


bool Polyhedron::globalIntersectionCheck()
{
    for (auto& edge : m_innerEdges)
        if (edgeGlobalIntersectionCheck(edge))
            return true;

    return false;
}




void Polyhedron::smoothMesh(unsigned iterationsNum)
{
    for (unsigned i = 0; i < iterationsNum; i++)
    {
        for (auto& vert : m_innerVerts)
        {
            Vec shift;
            int delta_shifts_num = 0;
            for (auto& edge : m_innerEdges)
            {
                if (vert == edge->verts[0])
                {
                    shift += *edge->verts[1] - *vert;
                    delta_shifts_num++;
                }
                else if (vert == edge->verts[1])
                {
                    shift += *edge->verts[0] - *vert;
                    delta_shifts_num++;
                }
            }
            shift /= delta_shifts_num;
            vert->setPos(vert->pos() + shift);
        }
    }
}


void Polyhedron::smoothNotFinisedMesh(unsigned iterationsNum)
{
    for (unsigned i = 0; i < iterationsNum; i++)
        for (auto &vert : m_innerVerts)
            smoothAroundFrontVert(vert);
}


void Polyhedron::smoothFront(unsigned iterationsNum)
{
    for (unsigned i = 0; i < iterationsNum; i++)
        for (auto &vert : m_innerVerts)
            smoothAroundFrontVert(vert);
}


void Polyhedron::smoothAroundFrontVert(Vert* fVert)
{
    if (fVert->belongsToShellFace ||
        fVert->belongsToShellEdge ||
        fVert->belongsToShellVert)
        return;

    Vec shift;
    int delta_shifts_num = 0;
    for (auto &edge : m_innerEdges)
    {
        Vec d_shift;
        if (fVert == edge->verts[0])
        {
            d_shift = *edge->verts[1] - *fVert;
        }
        else if (fVert == edge->verts[1])
        {
            d_shift = *edge->verts[0] - *fVert;
        }
        shift += d_shift * (d_shift.magnitude() - m_preferredLength);
        delta_shifts_num++;
    }
    shift /= delta_shifts_num;
    fVert->setPos(fVert->pos() + shift);
}


std::pair<real_t, real_t> Polyhedron::analyzeMeshQuality()
{
    size_t simps_num = 0;
    real_t av_q = 0.0;
    real_t min_q = 1.0;
    for (auto &tetr : m_innerTetrs)
    {
        real_t q = tetr->computeQuality();
        av_q += q;
        simps_num++;
        if (q < min_q)
            min_q = q;
    }
    av_q /= simps_num;

    return { min_q, av_q };
}




Polyhedron::Polyhedron() {}


Polyhedron::Polyhedron(PolyhedralSet* polyhedr)
{
    m_polycr = polyhedr;
}


Polyhedron::~Polyhedron()
{
    for (auto& tetr : m_innerTetrs)
        delete tetr;
    for (auto& face : m_innerFaces)
        delete face;
    for (auto& edge : m_innerEdges)
        delete edge;
    for (auto& vert : m_innerVerts)
        delete vert;
}




real_t Polyhedron::preferredLength() const
{
    return m_preferredLength;
}


void Polyhedron::addToShell(const shell::Face* shellFace)
{
    m_shellFaces.push_back(const_cast<shell::Face*>(shellFace));
}


void Polyhedron::addToShell(const shell::Edge* shellEdge)
{
    m_shellEdges.push_back(const_cast<shell::Edge*>(shellEdge));
}


void Polyhedron::addToShell(const shell::Vert* shellVert)
{
    m_shellVerts.push_back(const_cast<shell::Vert*>(shellVert));
}




bool Polyhedron::shellContains(const shell::Face* shellFace) const
{
    return std::find(m_shellFaces.begin(), m_shellFaces.end(), shellFace) != m_shellFaces.end();
}


bool Polyhedron::shellContains(const shell::Edge* shellEdge) const
{
    return std::find(m_shellEdges.begin(), m_shellEdges.end(), shellEdge) != m_shellEdges.end();
}


bool Polyhedron::shellContains(const shell::Vert* shellVert) const
{
    return std::find(m_shellVerts.begin(), m_shellVerts.end(), shellVert) != m_shellVerts.end();
}




const std::vector<Tetr*>& Polyhedron::innerTetrs() const
{
    return m_innerTetrs;
}


const std::vector<Face*>& Polyhedron::innerFaces() const
{
    return m_innerFaces;
}


const std::vector<Vert*>& Polyhedron::innerVerts() const
{
    return m_innerVerts;
}




const std::list<front::Face*>& Polyhedron::frontFaces() const
{
    return m_frontFaces;
}


const std::list<front::Edge*>& Polyhedron::frontEdges() const
{
    return m_frontEdges;
}
