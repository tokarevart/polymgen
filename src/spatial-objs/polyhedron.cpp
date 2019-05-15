// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/polyhedron.h"
#include <algorithm>
#include <iostream>
#include "helpers/spatial-algs/spatial-algs.h"


using namespace pmg;
using pair_rr = std::pair<real_t, real_t>;
using pair_ff = std::pair<front::Face*, front::Face*>;


#define DET(a, b, c, d) \
        (a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0[0], corner1[0], point[0]) && \
         BETWEEN(corner0[1], corner1[1], point[1]))


// TODO: reduce the amount of defines
#define ALPHA_P      static_cast<real_t>(70.52877936550931)
#define DEG_1_IN_RAD static_cast<real_t>( 0.0174532925199432957)

#define ONE_3                static_cast<real_t>(0.3333333333333333)
#define SQRT3_2              static_cast<real_t>(0.8660254037844386)
#define SQRT_2_3             static_cast<real_t>(0.8164965809277260)
#define ONE_PLUS_SQRT2_SQRT3 static_cast<real_t>(1.3938468501173517)

#define C_MIN_DIS              static_cast<real_t>(2e-1)
#define C_EDGES_INTERS_DIST static_cast<real_t>(4e-3)

#define C_MAXD static_cast<real_t>(0.3)
#define C_D    static_cast<real_t>(0.4)


template <typename T>
constexpr real_t degToRad( T value )
{
    return value * DEG_1_IN_RAD;
}


//#define DEV_DEBUG


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


bool Polyhedron::segmentIntersectMesh(const vec3& v0, const vec3& v1) const
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


bool Polyhedron::segmentIntersectFront(const vec3& v0, const vec3& v1) const
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


bool Polyhedron::edgeIntersectFront(const Vert* v0, const vec3& v1) const
{
    for (auto& fface : m_frontFaces)
    {
        if (relations::contains(fface, v0))
            continue;

        if (spatalgs::doesSegmentIntersectTriangle(
                v0->pos(), v1,
                fface->face->edges[0]->verts[0]->pos(),
                fface->face->edges[0]->verts[1]->pos(),
                fface->face->findVertNot(fface->face->edges[0])->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::edgeIntersectFront(const Vert* v0, const Vert* v1) const
{
    for (auto& fface : m_frontFaces)
    {
        if (   relations::contains(fface, v0)
            || relations::contains(fface, v1))
            continue;

        if (spatalgs::doesSegmentIntersectTriangle(
                v0->pos(), v1->pos(),
                fface->face->edges[0]->verts[0]->pos(),
                fface->face->edges[0]->verts[1]->pos(),
                fface->face->findVertNot(fface->face->edges[0])->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::edgeIntersectAnyFace(const Edge* edge) const
{
    for (auto& face : m_innerFaces)
    {
        if (   face->contains(edge->verts[0])
            || face->contains(edge->verts[1]))
            continue;

        if (spatalgs::doesSegmentIntersectTriangle(
                edge->verts[0]->pos(), edge->verts[1]->pos(),
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


bool Polyhedron::potentialEdgeIntersectFront(front::Edge* fEdge) const
{
    auto opp_verts = fEdge->oppVerts();

    if (edgeIntersectFront(std::get<0>(opp_verts), std::get<1>(opp_verts)))
        return true;

    for (auto& fedge : m_frontEdges)
        if (fedge->edge->contains(std::get<0>(opp_verts)) &&
            fedge->edge->contains(std::get<1>(opp_verts)))
            return false;

    for (auto& fedge : m_frontEdges)
    {
        Vert* vert_buf;
        std::array<bool, 2> contains =
        {
            fedge->edge->contains(std::get<0>(opp_verts)),
            fedge->edge->contains(std::get<1>(opp_verts))
        };
        // TODO: change that doubtful checks
        if (contains[0])
        {
            if (fedge->edge->verts[0] == std::get<0>(opp_verts))
                vert_buf = fedge->edge->verts[1];
            else
                vert_buf = fedge->edge->verts[0];

            if (spatalgs::distancePointToSegment(vert_buf->pos(), std::get<0>(opp_verts)->pos(), std::get<1>(opp_verts)->pos()) < C_EDGES_INTERS_DIST * m_prefLen)
                return true;
        }
        else if (contains[1])
        {
            if (fedge->edge->verts[0] == std::get<1>(opp_verts))
                vert_buf = fedge->edge->verts[1];
            else
                vert_buf = fedge->edge->verts[0];

            if (spatalgs::distancePointToSegment(vert_buf->pos(), std::get<0>(opp_verts)->pos(), std::get<1>(opp_verts)->pos()) < C_EDGES_INTERS_DIST * m_prefLen)
                return true;
        }
        else
        {
            if (spatalgs::segmentsDistance(
                    std::get<0>(opp_verts)->pos(), std::get<1>(opp_verts)->pos(),
                    fedge->edge->verts[0]->pos(), fedge->edge->verts[1]->pos()) < C_EDGES_INTERS_DIST * m_prefLen)
                return true;
        }
    }

    return false;
}


bool Polyhedron::anyEdgeIntersectFace(const Vert* v0, const Vert* v1, const vec3& v2) const
{
    for (auto& fedge : m_frontEdges)
    {
        std::array<bool, 2> contains =
        {
            fedge->edge->contains(v0),
            fedge->edge->contains(v1)
        };

        if (contains[0] || contains[1])
            continue;

        auto cur_verts = fedge->edge->verts;
        if (spatalgs::doesSegmentIntersectTriangle(
            cur_verts[0]->pos(), cur_verts[1]->pos(),
            v0->pos(), v1->pos(), v2))
            return true;
    }

    return false;
}


bool Polyhedron::anyEdgeIntersectFace(const Vert* v0, const Vert* v1, const Vert* v2) const
{
    for (auto& fedge : m_frontEdges)
    {
        if (   fedge->edge->contains(v0)
            || fedge->edge->contains(v1)
            || fedge->edge->contains(v2))
            continue;

        if (spatalgs::doesSegmentIntersectTriangle(
                fedge->edge->verts[0]->pos(), fedge->edge->verts[1]->pos(),
                v0->pos(), v1->pos(), v1->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::anyEdgeIntersectPotentialFaces(front::Edge* fEdge) const
{
    auto opp_verts = fEdge->oppVerts();
        
    if (anyEdgeIntersectFace(fEdge->edge->verts[0], std::get<0>(opp_verts), std::get<1>(opp_verts)) ||
        anyEdgeIntersectFace(fEdge->edge->verts[1], std::get<0>(opp_verts), std::get<1>(opp_verts)))
        return true;

    return false;
}


bool Polyhedron::anyVertInsidePotentialTetrCheck(front::Edge* fEdge) const
{
    auto opp_verts = fEdge->oppVerts();

    std::array<vec3, 4> points;
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
            spatalgs::isPointInTetrahedron(vert->pos(), points[0], points[1], points[2], points[3]))
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

    auto opp_fedge_opp_verts = opp_fedge->oppVerts();

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

    auto opp_ffaces = opp_fedge->adjFFaces();

    if (!fEdge->edge->contains(opp_ffaces.first->face->findVertNot(opp_fedge->edge)) ||
        !fEdge->edge->contains(opp_ffaces.second->face->findVertNot(opp_fedge->edge)))
        return false;

    return true;
}


bool Polyhedron::parallelFacesCheck(front::Edge* fEdge) const
{
    auto adj_ffaces = fEdge->adjFFaces();

    std::array<Vert*, 2> opp_verts;
    opp_verts[0] = std::get<0>(adj_ffaces)->face->findVertNot(fEdge->edge);
    opp_verts[1] = std::get<1>(adj_ffaces)->face->findVertNot(fEdge->edge);

    std::array<vec3, 2> plane0;
    plane0[0] = *opp_verts[0] - *fEdge->edge->verts[0];
    plane0[1] = *opp_verts[1] - *fEdge->edge->verts[0];
    vec3 normal0 = vec3::cross(plane0[0], plane0[1]).normalize();
    std::array<vec3, 2> plane1;
    plane1[0] = *opp_verts[0] - *fEdge->edge->verts[1];
    plane1[1] = *opp_verts[1] - *fEdge->edge->verts[1];
    vec3 normal1 = vec3::cross(plane1[0], plane1[1]).normalize();

    for (auto& fface : m_frontFaces)
    {
        std::array<Edge*, 2> inters_reses;
        if ((fface != std::get<0>(adj_ffaces)) &&
            (fface != std::get<1>(adj_ffaces)))
        {
            inters_reses[0] = relations::adjacentByEdge(fface->face, std::get<0>(adj_ffaces)->face);
            inters_reses[1] = relations::adjacentByEdge(fface->face, std::get<1>(adj_ffaces)->face);
            if (!XOR(static_cast<bool>(inters_reses[0]), static_cast<bool>(inters_reses[1])))
                continue;

            std::array<Vert*, 2> fface_to_verts;
            fface_to_verts[0] = fface->face->edges[0]->verts[0];
            fface_to_verts[1] = fface->face->edges[0]->verts[1];
            Vert* fface_from_vert = fface->face->findVertNot(fface->face->edges[0]);

            std::array<vec3, 2> f_plane;
            f_plane[0] = *fface_to_verts[0] - *fface_from_vert;
            f_plane[1] = *fface_to_verts[1] - *fface_from_vert;
            vec3 f_normal = vec3::cross(f_plane[0], f_plane[1]).normalize();

            if (std::abs(std::abs(vec3::dot(f_normal, normal0)) - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6) ||
                std::abs(std::abs(vec3::dot(f_normal, normal1)) - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6))
            {
                size_t i = inters_reses[0] ? 0 : 1;

                std::array<vec3, 2> border_verts;
                border_verts[0] = inters_reses[i]->verts[0]->pos();
                border_verts[1] = inters_reses[i]->verts[1]->pos();

                vec3 main_face_3rd_vert;
                if (inters_reses[i]->contains(opp_verts[0]))
                    main_face_3rd_vert = opp_verts[1]->pos();
                else
                    main_face_3rd_vert = opp_verts[0]->pos();

                vec3 cur_face_3rd_vert = fface->face->findVertNot(inters_reses[i])->pos();

                vec3 main_face_cross = vec3::cross(main_face_3rd_vert - border_verts[0], main_face_3rd_vert - border_verts[1]);
                vec3 cur_face_cross  = vec3::cross(cur_face_3rd_vert  - border_verts[0], cur_face_3rd_vert  - border_verts[1]);
                if (vec3::dot(main_face_cross, cur_face_cross) > static_cast<real_t>(0.0))
                    return true;
            }
        }
    }

    return false;
}


bool Polyhedron::doesFrontIntersectSphere(const vec3& center, real_t radius) const
{
    for (auto& fface : m_frontFaces)
    {
        std::array<vec3, 3> triangle;
        triangle[0] = fface->face->edges[0]->verts[0]->pos();
        triangle[1] = fface->face->edges[0]->verts[1]->pos();
        triangle[2] = fface->face->findVertNot(fface->face->edges[0])->pos();
        if (spatalgs::doesTriangleIntersectSphere(triangle[0], triangle[1], triangle[2], center, radius))
            return true;
    }

    return false;
}


std::pair<real_t, real_t> Polyhedron::computeMinMaxEdgesLengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3)
{
    auto min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    min_max.first  = std::sqrt(min_max.first);
    min_max.second = std::sqrt(min_max.second);
    return min_max;
}


std::pair<real_t, real_t> Polyhedron::computeMinMaxEdgesSqrLengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3)
{
    std::array<real_t, 6> sqr_magns;
    sqr_magns[0] = (p1 - p0).sqrMagnitude();
    sqr_magns[1] = (p2 - p0).sqrMagnitude();
    sqr_magns[2] = (p3 - p0).sqrMagnitude();
    sqr_magns[3] = (p2 - p1).sqrMagnitude();
    sqr_magns[4] = (p3 - p1).sqrMagnitude();
    sqr_magns[5] = (p3 - p2).sqrMagnitude();
    return std::minmax({ sqr_magns[0], sqr_magns[1], sqr_magns[2], sqr_magns[3], sqr_magns[4], sqr_magns[5] });
}


real_t Polyhedron::computeTetrSimpleQuality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return std::sqrt(sqr_min_max.first / sqr_min_max.second);
}


real_t Polyhedron::computeTetrSimpleSqrQuality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3)
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
    if (curFEdge->angle() < degToRad(70))
        return true;

    auto opp_verts = curFEdge->oppVerts();
    if (curFEdge->findOppEdge() ||
        (curFEdge->angle() > degToRad(70) &&
         curFEdge->angle() < degToRad(100) &&
         (*std::get<1>(opp_verts) - *std::get<0>(opp_verts)).sqrMagnitude() <= m_prefLen * m_prefLen))
        return true;

    if (frontSplitCheck(curFEdge))
        return true;

    return false;
}


bool Polyhedron::exhaustWithNewVertPriorityPredicate(front::Edge* currentFrontEdge)
{
    if (currentFrontEdge->angle() > degToRad(120))
        return true;

    return false;
}

 
Polyhedron::ExhaustType Polyhedron::computeExhaustionTypeQualityPriority(
    front::Edge* currentFrontEdge,
    front::Face*& out_withNWFrontFace, vec3*& out_withNWNewVertPos)
{
    if (frontSplitCheck(currentFrontEdge))
        return ExhaustType::WithoutNewVert;
    
    if (parallelFacesCheck(currentFrontEdge) ||
        potentialEdgeIntersectFront(currentFrontEdge) ||
        anyEdgeIntersectPotentialFaces(currentFrontEdge) ||
        anyVertInsidePotentialTetrCheck(currentFrontEdge))
        return ExhaustType::WithNewVert;

    auto opp_verts = currentFrontEdge->oppVerts();
    real_t without_nv_quality = computeTetrSimpleSqrQuality(
        currentFrontEdge->edge->verts[0]->pos(),
        currentFrontEdge->edge->verts[1]->pos(),
        std::get<0>(opp_verts)->pos(),
        std::get<1>(opp_verts)->pos());

    front::Face* fface = chooseFaceForExhaustionWithNewVert(currentFrontEdge);
    vec3 new_vert_pos;
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
    out_withNWNewVertPos = new vec3(new_vert_pos);
    return ExhaustType::WithNewVert;
}


vec3 Polyhedron::computeNormalInTetr(const front::Face* fFace, const vec3& oppVertPos) const
{
    return computeNormalInTetr(fFace->face->edges[0]->verts[0]->pos(),
                               fFace->face->edges[0]->verts[1]->pos(),
                               fFace->face->findVertNot(fFace->face->edges[0])->pos(),
                               oppVertPos);
}


vec3 Polyhedron::computeNormalInTetr(const front::Face* fFace, const pmg::Edge* oneOfRemainingEdges) const
{
    vec3 opp_pos = fFace->face->contains(oneOfRemainingEdges->verts[0]) ?
        oneOfRemainingEdges->verts[1]->pos() :
        oneOfRemainingEdges->verts[0]->pos();

    return computeNormalInTetr(fFace, opp_pos);
}


vec3 Polyhedron::computeNormalInTetr(const vec3& fFacePos0, const vec3& fFacePos1, const vec3& fFacePos2, const vec3& oppVertPos) const
{
    vec3 normal = vec3::cross(
        fFacePos0 - fFacePos2,
        fFacePos1 - fFacePos2).normalize();
    if (vec3::dot(normal, oppVertPos - fFacePos2) > static_cast<real_t>(0.0))
        normal *= static_cast<real_t>(-1.0);

    return normal;
}


void Polyhedron::setFEdgesInFrontSplit(const front::Edge* fEdge, std::array<front::Edge*, 2> newOppFEdges, std::array<front::Face*, 2> newFFaces, pair_ff oppFFaces) const
{
    pmg::Edge* opp_edge = newOppFEdges[0]->edge;
    std::array<vec3, 2> opp_verts_poses;
    opp_verts_poses[0] = opp_edge->verts[0]->pos();
    opp_verts_poses[1] = opp_edge->verts[1]->pos();
    vec3 main_vert0_pos  = newFFaces[0]->face->contains(fEdge->edge->verts[0]) ?
                fEdge->edge->verts[0]->pos() : fEdge->edge->verts[1]->pos();
    vec3 main_vert1_pos  = newFFaces[1]->face->contains(fEdge->edge->verts[0]) ?
                fEdge->edge->verts[0]->pos() : fEdge->edge->verts[1]->pos();
    vec3 main_vert0_proj = spatalgs::project(main_vert0_pos, opp_verts_poses[0], opp_verts_poses[1]);
    vec3 main_vert1_proj = spatalgs::project(main_vert1_pos, opp_verts_poses[0], opp_verts_poses[1]);
    vec3 main_vec0 = main_vert0_pos - main_vert0_proj;
    vec3 main_vec1 = main_vert1_pos - main_vert1_proj;

    vec3 adj_opp_pos0 = oppFFaces.first->face->findVertNot(opp_edge)->pos();
    vec3 adj_opp_pos1 = oppFFaces.second->face->findVertNot(opp_edge)->pos();
    vec3 adj_opp_proj0 = spatalgs::project(adj_opp_pos0, opp_verts_poses[0], opp_verts_poses[1]);
    vec3 adj_opp_proj1 = spatalgs::project(adj_opp_pos1, opp_verts_poses[0], opp_verts_poses[1]);
    vec3 adj_vec0 = adj_opp_pos0 - adj_opp_proj0;
    vec3 adj_vec1 = adj_opp_pos1 - adj_opp_proj1;

    std::array<std::array<real_t, 2>, 2> coses;
    coses[0][0] = vec3::cos(main_vec0, adj_vec0);
    coses[0][1] = vec3::cos(main_vec0, adj_vec1);
    coses[1][0] = vec3::cos(main_vec1, adj_vec0);
    coses[1][1] = vec3::cos(main_vec1, adj_vec1);

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

    size_t best_ofi;
    if (coses[0][0] > coses[0][1] && coses[1][0] > coses[1][1])
        best_ofi = 0;
    else
        // When coses[0][1] >= coses[0][0] && coses[1][1] >= coses[1][0]
        best_ofi = 1;

    size_t worst_ofi = best_ofi == 0 ? 1 : 0;

    std::array<front::Face*, 2> opp_ffaces = { oppFFaces.first, oppFFaces.second };
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


void Polyhedron::exhaustFrontCollapse(front::Edge* fEdge, front::Edge* oppFEdge)
{
    auto fedge_adj_faces = fEdge->adjFFaces();
    auto opp_ffaces   = oppFEdge->adjFFaces();

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
//        m_polyhset->output(PolyhedralSet::FileType::WavefrontObj, "debug.obj");


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
    // HACK: this volatile helps to avoid computational error
    volatile auto adj_ffaces = fEdge->adjFFaces();

    std::array<Vert*, 2> opp_verts;
    opp_verts[0] = adj_ffaces.first->face->findVertNot(fEdge->edge);
    opp_verts[1] = adj_ffaces.second->face->findVertNot(fEdge->edge);

    auto opp_edge   = oppFEdge->edge;
    auto opp_ffaces = oppFEdge->adjFFaces();

    opp_ffaces.first->removeFEdge(oppFEdge);
    opp_ffaces.second->removeFEdge(oppFEdge);
    removeFromFront(oppFEdge);
    std::array<front::Edge*, 2> new_opp_fedges;
    new_opp_fedges[0] = addToFront(opp_edge, false);
    new_opp_fedges[1] = addToFront(opp_edge, false);

    std::array<front::Edge*, 3> new_tetr_fedges;
    new_tetr_fedges[0] = nullptr;
    new_tetr_fedges[1] = nullptr;
    new_tetr_fedges[2] = new_opp_fedges[0];

    std::array<front::Face*, 2> new_ffaces;

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


    setFEdgesInFrontSplit(fEdge, { new_opp_fedges[0], new_opp_fedges[1] }, { new_ffaces[0], new_ffaces[1] }, opp_ffaces);

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
//        m_polyhset->output(PolyhedralSet::FileType::WavefrontObj, "debug.obj");


    removeFromFront(adj_ffaces.first);
    removeFromFront(adj_ffaces.second);
    removeFromFront(fEdge);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustFrontSplit";
#endif
}


void Polyhedron::exhaustWithoutNewVertOppEdgeExists(front::Edge* fEdge, front::Edge* oppFEdge)
{
    auto opp_ffaces = oppFEdge->adjFFaces();

    std::array<front::Face*, 3> main_ffaces;
    Vert* main_vert;
    auto fedge_adj_faces = fEdge->adjFFaces();
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

    std::array<front::Edge*, 3> new_tetr_fedges;
    for (size_t i = 0; i < 3; i++)
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
//        m_polyhset->output(PolyhedralSet::FileType::WavefrontObj, "debug.obj");

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
    // HACK: This volatile helps to avoid computational error.
    volatile auto adj_ffaces = fEdge->adjFFaces();

    std::array<Vert*, 2> opp_verts;
    opp_verts[0] = adj_ffaces.first->face->findVertNot(fEdge->edge);
    opp_verts[1] = adj_ffaces.second->face->findVertNot(fEdge->edge);
    
    front::Edge* opp_fedge = addToFront(new Edge(opp_verts[0], opp_verts[1]));

    std::array<front::Edge*, 3> new_tetr_fedges;
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
//        m_polyhset->output(PolyhedralSet::FileType::WavefrontObj, "debug.obj");

    
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
            exhaustFrontCollapse(fEdge, opp_fedge);

        else if (frontSplitCheck(fEdge, opp_fedge))
            exhaustFrontSplit(fEdge, opp_fedge);

        else
            exhaustWithoutNewVertOppEdgeExists(fEdge, opp_fedge);
    }
    else
    {
        exhaustWithoutNewVertOppEdgeDontExists(fEdge);
    }
}


bool Polyhedron::tryComputeNewVertPosType3(front::Face* fFace, vec3& out_pos)
{
    std::array<front::Edge*, 2> main_fedges;
    main_fedges[0] = fFace->fEdges[0];
    main_fedges[1] = fFace->fEdges[1];
    std::array<Edge*, 2> main_edges;
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

    auto adj_ffaces = main_fedges[0]->adjFFaces();
    auto fn0 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    adj_ffaces = main_fedges[1]->adjFFaces();
    auto fn1 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    adj_ffaces = third_f_edge->adjFFaces();
    auto fn2 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);

    vec3 v0_pos = v0->pos();
    vec3 v1_pos = v1->pos();
    vec3 v2_pos = v2->pos();

    vec3 n_m  = fFace->normal;
    vec3 n_n0 = fn0->normal;
    vec3 n_n1 = fn1->normal;
    vec3 n_n2 = fn2->normal;

    vec3 e_mn0 = (n_m + n_n0).normalize();
    vec3 e_mn1 = (n_m + n_n1).normalize();
    vec3 e_mn2 = n_m + n_n2;

    vec3 np_mn0 = vec3::cross(v2_pos - v0_pos, e_mn0);
    vec3 np_mn1 = vec3::cross(v2_pos - v1_pos, e_mn1);

    vec3 e = vec3::cross(np_mn0, np_mn1);

    vec3 new_pos = spatalgs::lineIntersectPlane(v2_pos, e, v0_pos, v1_pos, v0_pos + e_mn2);

    std::array<real_t, 4> sum_magns;
    size_t i = 0;
    for (auto& fface : { fn0, fn1, fn2, fFace })
        sum_magns[i++] =  fface->face->edges[0]->magnitude()
                        + fface->face->edges[1]->magnitude()
                        + fface->face->edges[2]->magnitude();
    real_t av_magn = (sum_magns[0] + sum_magns[1] + sum_magns[2] + sum_magns[3]) / static_cast<real_t>(12.0);
    if (edgeIntersectFront(v0, new_pos) ||
        edgeIntersectFront(v1, new_pos) ||
        edgeIntersectFront(v2, new_pos) ||
        doesFrontIntersectSphere(new_pos, C_MIN_DIS * av_magn) ||
        anyEdgeIntersectFace(v0, v1, new_pos) ||
        anyEdgeIntersectFace(v0, v2, new_pos) ||
        anyEdgeIntersectFace(v1, v2, new_pos))
        return false;

//    std::cout << std::endl << "Type3";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::tryComputeNewVertPosType2(front::Face* fFace, vec3& out_pos, size_t smallAngleIdx0, size_t smallAngleIdx1)
{
    std::array<front::Edge*, 2> main_fedges;
    main_fedges[0] = fFace->fEdges[smallAngleIdx0];
    main_fedges[1] = fFace->fEdges[smallAngleIdx1];
    std::array<Edge*, 2> main_edges;
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

    auto adj_ffaces = main_fedges[0]->adjFFaces();
    auto fn0 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn0 = fn0->face->findVertNot(main_fedges[0]->edge);
    adj_ffaces = main_fedges[1]->adjFFaces();
    auto fn1 = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn1 = fn1->face->findVertNot(main_fedges[1]->edge);

    vec3 v0_pos = v0->pos();
    vec3 v1_pos = v1->pos();
    vec3 v2_pos = v2->pos();
    vec3 vn0_pos = vn0->pos();
    vec3 vn1_pos = vn1->pos();

    vec3 n_m = fFace->normal;
    vec3 n_n0 = fn0->normal;
    vec3 n_n1 = fn1->normal;
    
    vec3 e_mn0 = (n_m + n_n0).normalize();
    vec3 e_mn1 = (n_m + n_n1).normalize();

    vec3 np_mn0 = vec3::cross(v2_pos - v0_pos, e_mn0);
    vec3 np_mn1 = vec3::cross(v2_pos - v1_pos, e_mn1);

    vec3 e = vec3::cross(np_mn0, np_mn1).normalize();
    if (vec3::dot(e, n_m) < static_cast<real_t>(0.0)) e *= static_cast<real_t>(-1.0);

    real_t lm0 = main_edges[0]->magnitude();
    real_t lm1 = main_edges[1]->magnitude();
    real_t l0 = (v1_pos - v0_pos).magnitude();
    real_t l1 = (vn0_pos - v0_pos).magnitude();
    real_t l2 = (vn0_pos - v2_pos).magnitude();
    real_t l3 = (vn1_pos - v1_pos).magnitude();
    real_t l4 = (vn1_pos - v2_pos).magnitude();
    real_t av_magn = (lm0 + lm1 + l0 + l1 + l2 + l3 + l4) / static_cast<real_t>(7.0);
    real_t raw_deform = C_D * (m_prefLen - av_magn);
    real_t deform = raw_deform < av_magn * C_MAXD ? raw_deform : av_magn * C_MAXD;
    real_t magn_d = av_magn + deform;
    vec3 new_pos = v2_pos + magn_d * e;

    if (edgeIntersectFront(v0, new_pos) ||
        edgeIntersectFront(v1, new_pos) ||
        edgeIntersectFront(v2, new_pos) ||
        doesFrontIntersectSphere(new_pos, C_MIN_DIS * magn_d) ||
        anyEdgeIntersectFace(v0, v1, new_pos) ||
        anyEdgeIntersectFace(v0, v2, new_pos) ||
        anyEdgeIntersectFace(v1, v2, new_pos))
    {
        // NOTE: do i really need it?

        new_pos = new_pos = v2_pos + av_magn * e;
        if (edgeIntersectFront(v0, new_pos) ||
            edgeIntersectFront(v1, new_pos) ||
            edgeIntersectFront(v2, new_pos) ||
            doesFrontIntersectSphere(new_pos, C_MIN_DIS * av_magn) ||
            anyEdgeIntersectFace(v0, v1, new_pos) ||
            anyEdgeIntersectFace(v0, v2, new_pos) ||
            anyEdgeIntersectFace(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type2";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::tryComputeNewVertPosType1(front::Face* fFace, vec3& out_pos, size_t smallAngleIdx)
{
    auto main_f_edge = fFace->fEdges[smallAngleIdx];
    auto main_edge   = fFace->fEdges[smallAngleIdx]->edge;
    std::array<vec3, 2> main_edge_poses;
    main_edge_poses[0] = main_edge->verts[0]->pos();
    main_edge_poses[1] = main_edge->verts[1]->pos();

    auto v0 = main_edge->verts[0];
    auto v1 = main_edge->verts[1];
    auto v2 = fFace->face->findVertNot(main_edge);

    auto adj_ffaces = main_f_edge->adjFFaces();
    auto fn = std::get<0>(adj_ffaces) == fFace ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn = fn->face->findVertNot(main_edge);

    vec3 v0_pos = v0->pos();
    vec3 v1_pos = v1->pos();
    vec3 v2_pos = v2->pos();
    vec3 vn_pos = vn->pos();

    vec3 v2pr = spatalgs::project(v2_pos, v0_pos, v1_pos);
    vec3 vnpr = spatalgs::project(vn_pos, v0_pos, v1_pos);

    vec3 c = static_cast<real_t>(0.25) * (v0_pos + v1_pos + v2pr + vnpr);
    vec3 e = (fFace->normal + fn->normal).normalize();
    real_t me_magn = main_edge->magnitude();
    real_t l0 = (v2_pos - v0_pos).magnitude();
    real_t l1 = (v2_pos - v1_pos).magnitude();
    real_t l2 = (vn_pos - v0_pos).magnitude();
    real_t l3 = (vn_pos - v1_pos).magnitude();
    real_t av_magn = static_cast<real_t>(0.2) * (me_magn + l0 + l1 + l2 + l3);
    real_t v0_c_dist = (c - v0_pos).magnitude();
    real_t raw_deform = C_D * (m_prefLen - av_magn);
    real_t deform = raw_deform < av_magn * C_MAXD ? raw_deform : av_magn * C_MAXD;
    real_t magn_d = av_magn + deform;
    vec3 new_pos = c + std::sqrt(magn_d * magn_d - v0_c_dist * v0_c_dist) * e;

    if (edgeIntersectFront(v0, new_pos) ||
        edgeIntersectFront(v1, new_pos) ||
        edgeIntersectFront(v2, new_pos) ||
        doesFrontIntersectSphere(new_pos, C_MIN_DIS * magn_d) ||
        anyEdgeIntersectFace(v0, v1, new_pos) ||
        anyEdgeIntersectFace(v0, v2, new_pos) ||
        anyEdgeIntersectFace(v1, v2, new_pos))
    {
        // NOTE: do i really need it?

        new_pos = c + std::sqrt(av_magn * av_magn - v0_c_dist * v0_c_dist) * e;
        if (edgeIntersectFront(v0, new_pos) ||
            edgeIntersectFront(v1, new_pos) ||
            edgeIntersectFront(v2, new_pos) ||
            doesFrontIntersectSphere(new_pos, C_MIN_DIS * av_magn) ||
            anyEdgeIntersectFace(v0, v1, new_pos) ||
            anyEdgeIntersectFace(v0, v2, new_pos) ||
            anyEdgeIntersectFace(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type1";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::tryComputeNewVertPosType0(front::Face* fFace, vec3& out_pos)
{
    real_t av_magn = ONE_3 * (  fFace->face->edges[0]->magnitude()
                              + fFace->face->edges[1]->magnitude()
                              + fFace->face->edges[2]->magnitude());
    real_t raw_deform = C_D * (m_prefLen - av_magn);
    real_t deform = raw_deform < av_magn * C_MAXD ? raw_deform : av_magn * C_MAXD;
    real_t magn_d = av_magn + deform;
    vec3 new_pos = fFace->computeCenter() + std::sqrt(magn_d * magn_d - ONE_3 * av_magn * av_magn) * fFace->normal;

    auto v0 = fFace->face->edges[0]->verts[0];
    auto v1 = fFace->face->edges[0]->verts[1];
    auto v2 = fFace->face->findVertNot(fFace->face->edges[0]);

    if (edgeIntersectFront(v0, new_pos) ||
        edgeIntersectFront(v1, new_pos) ||
        edgeIntersectFront(v2, new_pos) ||
        doesFrontIntersectSphere(new_pos, C_MIN_DIS * magn_d) ||
        anyEdgeIntersectFace(v0, v1, new_pos) ||
        anyEdgeIntersectFace(v0, v2, new_pos) ||
        anyEdgeIntersectFace(v1, v2, new_pos))
    {
        // NOTE: do i really need it?

        new_pos = fFace->computeCenter() + std::sqrt(av_magn * av_magn - ONE_3 * av_magn * av_magn) * fFace->normal;
        if (edgeIntersectFront(v0, new_pos) ||
            edgeIntersectFront(v1, new_pos) ||
            edgeIntersectFront(v2, new_pos) ||
            doesFrontIntersectSphere(new_pos, C_MIN_DIS * av_magn) ||
            anyEdgeIntersectFace(v0, v1, new_pos) ||
            anyEdgeIntersectFace(v0, v2, new_pos) ||
            anyEdgeIntersectFace(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type0";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::tryComputeNewVertPos(front::Face* fFace, vec3& out_pos)
{
    std::array<real_t, 3> angs =
    { 
        fFace->fEdges[0]->angle(),
        fFace->fEdges[1]->angle(),
        fFace->fEdges[2]->angle()
    };
    std::array<size_t, 3> idces;
    size_t n_small_angs = 0;
    if (angs[0] < degToRad(140)) idces[n_small_angs++] = 0;
    if (angs[1] < degToRad(140)) idces[n_small_angs++] = 1;
    if (angs[2] < degToRad(140)) idces[n_small_angs++] = 2;

    switch (n_small_angs)
    {
    case 0: return tryComputeNewVertPosType0(fFace, out_pos);
    case 1: return tryComputeNewVertPosType1(fFace, out_pos, idces[0]);
    case 2: return tryComputeNewVertPosType2(fFace, out_pos, idces[0], idces[1]);
    case 3: return tryComputeNewVertPosType3(fFace, out_pos);
    }

    return true;
}


real_t Polyhedron::sqr4FaceArea(const front::Face* fFace) const
{
    vec3 edge0_vec = *fFace->face->edges[0]->verts[1] - *fFace->face->edges[0]->verts[0];
    vec3 edge1_vec = *fFace->face->edges[1]->verts[1] - *fFace->face->edges[1]->verts[0];
    return vec3::cross(edge0_vec, edge1_vec).sqrMagnitude();
}


front::Face* Polyhedron::chooseFaceForExhaustionWithNewVert(front::Edge* fEdge)
{
    auto adj_ffaces = fEdge->adjFFaces();

    return sqr4FaceArea(std::get<0>(adj_ffaces)) < sqr4FaceArea(std::get<1>(adj_ffaces)) ?
//    return std::get<0>(adj_ffaces)->computeQuality() < std::get<1>(adj_ffaces)->computeQuality() ?
        std::get<0>(adj_ffaces) :
        std::get<1>(adj_ffaces);
}


void Polyhedron::exhaustWithNewVert(front::Face* fFace, const vec3& vertPos)
{
    Vert* new_vert = new Vert(vertPos);
    m_innerVerts.push_back(new_vert);

    std::array<front::Edge*, 6> new_tetr_fedges;
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
//        m_polyhset->output(PolyhedralSet::FileType::WavefrontObj, "debug.obj");

    removeFromFront(fFace);
}


bool Polyhedron::tryExhaustWithoutNewVert(front::Edge* fEdge)
{
    // TODO: improve that checks
    if (parallelFacesCheck(fEdge) ||
        anyEdgeIntersectPotentialFaces(fEdge) ||
        anyVertInsidePotentialTetrCheck(fEdge))
        return false;

    exhaustWithoutNewVert(fEdge);
    return true;
}


bool Polyhedron::tryExhaustWithNewVert(front::Edge* frontEdge)
{
    if (parallelFacesCheck(frontEdge))
        return false;

    auto exhaust_fface = chooseFaceForExhaustionWithNewVert(frontEdge);
    vec3 new_vert_pos;
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
            vec3*        new_vert_pos       = nullptr;
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
//        m_polyhset->output(PolyhedralSet::FileType::WavefrontObj, "debug.obj");
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
    vec3 accum;
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

    std::cout << "\n{ " << accum.x[0] + accum.x[1] + accum.x[2] << " }";
}


void Polyhedron::tetrahedralize(real_t preferredLen, genparams::Volume genParams)
{
    m_prefLen = preferredLen;
    initializeFront();
    computeFrontNormals();
    processAngles();
//    if (globalIntersectionCheck())
//        throw std::logic_error("Intersection error.\npmg::Polyhedron::globalIntersectionCheck returned true.");

    smoothMesh(genParams.nSmoothIters);

    m_isQualityAnalyzed     = false;
}


bool Polyhedron::globalIntersectionCheck()
{
    for (auto& edge : m_innerEdges)
        if (edgeIntersectAnyFace(edge))
            return true;

    return false;
}


void Polyhedron::smoothMesh(size_t nIters)
{
    for (size_t i = 0; i < nIters; i++)
    {
        for (auto& vert : m_innerVerts)
        {
            vec3 shift;
            size_t n_delta_shifts = 0;
            for (auto& edge : m_innerEdges)
            {
                if (vert == edge->verts[0])
                {
                    shift += *edge->verts[1] - *vert;
                    n_delta_shifts++;
                }
                else if (vert == edge->verts[1])
                {
                    shift += *edge->verts[0] - *vert;
                    n_delta_shifts++;
                }
            }
            shift /= n_delta_shifts;
            vert->setPos(vert->pos() + shift);
        }
    }
}


void Polyhedron::smoothNotFinisedMesh(size_t nIters)
{
    for (size_t i = 0; i < nIters; i++)
        for (auto &vert : m_innerVerts)
            smoothAroundFrontVert(vert);
}


void Polyhedron::smoothFront(size_t nIters)
{
    for (size_t i = 0; i < nIters; i++)
        for (auto &vert : m_innerVerts)
            smoothAroundFrontVert(vert);
}


void Polyhedron::smoothAroundFrontVert(Vert* fVert)
{
    if (fVert->belongsToSFace ||
        fVert->belongsToSEdge ||
        fVert->belongsToSVert)
        return;

    vec3 shift;
    size_t n_delta_shifts = 0;
    for (auto &edge : m_innerEdges)
    {
        vec3 d_shift;
        if (fVert == edge->verts[0])
        {
            d_shift = *edge->verts[1] - *fVert;
        }
        else if (fVert == edge->verts[1])
        {
            d_shift = *edge->verts[0] - *fVert;
        }
        shift += d_shift * (d_shift.magnitude() - m_prefLen);
        n_delta_shifts++;
    }
    shift /= n_delta_shifts;
    fVert->setPos(fVert->pos() + shift);
}


pair_rr Polyhedron::analyzeMeshQuality(std::list<Tetr*>::iterator* out_minQualityTetr)
{
    if (m_isQualityAnalyzed)
        return m_meshQuality;

    size_t n_tetrs = 0;
    real_t av_q = static_cast<real_t>(0.0);
    real_t min_q = static_cast<real_t>(1.0);
    std::list<Tetr*>::iterator min_q_tetr = m_innerTetrs.begin();
    for (auto iter = m_innerTetrs.begin(); iter != m_innerTetrs.end(); iter++)
    {
        real_t q = (*iter)->computeQuality();
        av_q += q;
        n_tetrs++;
        if (q < min_q)
        {
            min_q = q;
            if (out_minQualityTetr)
                min_q_tetr = iter;
        }
    }
    av_q /= n_tetrs;
    if (out_minQualityTetr)
        *out_minQualityTetr = min_q_tetr;

    m_isQualityAnalyzed = true;
    return m_meshQuality = { min_q, av_q };
}


pair_rr Polyhedron::analyzeMeshAbsGrad()
{
    if (m_isMeshAbsGradAnalyzed)
        return m_meshAbsGrad;

    for (auto& tetr : m_innerTetrs)
    {
        real_t volume = tetr->computeVolume();
        for (auto& vert : tetr->verts)
        {
            vert->minAdjTetrVol = std::min(volume, vert->minAdjTetrVol);
            vert->maxAdjTetrVol = std::max(volume, vert->maxAdjTetrVol);
        }
    }

    real_t min_abs_grad = static_cast<real_t>(1.0);
    real_t av_abs_grad  = static_cast<real_t>(0.0);
    for (auto& vert : m_innerVerts)
    {
        real_t cur_abs_grad = vert->minAdjTetrVol / vert->maxAdjTetrVol;
        min_abs_grad = std::min(cur_abs_grad, min_abs_grad);
        av_abs_grad += cur_abs_grad;
    }
    av_abs_grad /= m_innerVerts.size();

    m_isMeshAbsGradAnalyzed = true;
    return m_meshAbsGrad = { min_abs_grad, av_abs_grad };
}




Polyhedron::Polyhedron() {}


Polyhedron::Polyhedron(PolyhedralSet* polyhset)
{
    m_polyhset = polyhset;
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
    return m_prefLen;
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


const std::list<Tetr*>& Polyhedron::innerTetrs() const
{
    return m_innerTetrs;
}


const std::list<Face*>& Polyhedron::innerFaces() const
{
    return m_innerFaces;
}


const std::list<Vert*>& Polyhedron::innerVerts() const
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
