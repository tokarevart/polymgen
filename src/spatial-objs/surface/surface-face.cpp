// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/surface/surface-face.h"
#include <algorithm>
#include "helpers/spatial-algs/spatial-algs.h"


using namespace pmg;
namespace sfront = surface::front;
using pair_rr = std::pair<real_t, real_t>;
using pair_ff = std::pair<pmg::Face*, pmg::Face*>;
using pair_ee = std::pair<pmg::Edge*, pmg::Edge*>;


#define DEG_1_IN_RAD static_cast<real_t>(0.0174532925199432957)
#define PI           static_cast<real_t>(M_PI)

#define NINE_DIV_SIXTEEN static_cast<real_t>(0.5625)
#define SIXTEEN_DIV_NINE static_cast<real_t>(1.7777777777777777)
#define SQRT3_2          static_cast<real_t>(0.8660254037844386)

#define K_MIN_DIS              static_cast<real_t>(2e-1)
#define C_EDGES_INTERS_DIST static_cast<real_t>(1e-4)

#define K_MAXD static_cast<real_t>(0.3)
#define K_D    static_cast<real_t>(0.5)


template <typename T>
constexpr real_t degToRad(T value)
{
    return value * DEG_1_IN_RAD;
}


real_t surface::Face::preferredLength() const
{
    return m_prefLen;
}


const std::list<pmg::Face*>& surface::Face::innerFaces() const
{
    return m_innerFaces;
}


const std::list<pmg::Edge*>& surface::Face::innerEdges() const
{
    return m_innerEdges;
}


const std::list<pmg::Vert*>& surface::Face::innerVerts() const
{
    return m_innerVerts;
}


const std::list<sfront::Edge*>& surface::Face::frontEdges() const
{
    return m_frontEdges;
}


surface::Vert* surface::Face::findVertNot(const surface::Edge* edge) const
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


surface::Edge* surface::Face::findSurfaceEdgeContaining(const pmg::Edge* edge) const
{
    for (auto& sedge : edges)
        if (sedge->contains(edge))
            return sedge;

    return nullptr;
}


void surface::Face::triangulate(real_t preferredLen, genparams::Surface genParams)
{
    m_prefLen = preferredLen;
    initializeFront();
    computeFrontNormals();
    processAngles();
//    if (globalIntersectionCheck())
//        throw std::logic_error("Intersection error.\npmg::surface::Face::globalIntersectionCheck returned true.");

    optimizeMesh(genParams.nSmoothIters, genParams.nDelaunaySmoothIters);
}


bool surface::Face::contains(const surface::Edge* edge) const
{
    for (auto& edge_ : edges)
        if (edge_ == edge)
            return true;

    return false;
}


bool surface::Face::contains(const surface::Vert* vert) const
{
    for (auto& edge : edges)
        if (edge->contains(vert))
            return true;

    return false;
}




surface::Face::Face(
    const surface::Edge* edge0,
    const surface::Edge* edge1,
    const surface::Edge* edge2)
{
    edges[0] = const_cast<surface::Edge*>(edge0);
    edges[1] = const_cast<surface::Edge*>(edge1);
    edges[2] = const_cast<surface::Edge*>(edge2);
}




sfront::Vert* surface::Face::findFrontVert(const pmg::Vert* vert) const
{
    for (auto& fvert : m_frontVerts)
        if (fvert->vert == vert)
            return fvert;

    return nullptr;
}


sfront::Edge* surface::Face::addToFront(const pmg::Edge* edge)
{
    sfront::Edge* new_f_edge = new sfront::Edge(this, edge);
    m_frontEdges.push_back(new_f_edge);
    m_innerEdges.push_back(new_f_edge->edge);
    return new_f_edge;
}


sfront::Vert* surface::Face::addToFront(const pmg::Vert* vert)
{
    sfront::Vert* new_f_vert = new sfront::Vert(this, vert);
    m_frontVerts.push_back(new_f_vert);
    m_innerVerts.push_back(new_f_vert->vert);
    return new_f_vert;
}


void surface::Face::removeFromFront(sfront::Edge* fEdge)
{
    m_frontEdges.erase(std::find(m_frontEdges.begin(), m_frontEdges.end(), fEdge));
    delete fEdge;
}


void surface::Face::removeFromFront(sfront::Vert* fVert)
{
    m_frontVerts.erase(std::find(m_frontVerts.begin(), m_frontVerts.end(), fVert));
    delete fVert;
}


bool surface::Face::anyVertInsidePotentialTriangCheck(sfront::Vert* fVert) const
{
    auto opp_verts = fVert->oppVerts();
    std::array<pmg::Vert*, 3> tr
    {
        fVert->vert,
        opp_verts.first,
        opp_verts.second
    };

    for (auto& fvert : m_frontVerts)
        if (fvert->vert != tr[0] && fvert->vert != tr[1] && fvert->vert != tr[2] &&
            spatalgs::isPointOnTriangle(fvert->vert->pos(), tr[0]->pos(), tr[1]->pos(), tr[2]->pos()))
            return true;

    return false;
}


bool surface::Face::doesSegmentIntersectsWithFront(const vec3& v0, const vec3& v1) const
{
    for (auto& fedge : m_frontEdges)
        if (spatalgs::segmentsDistance(
                fedge->edge->verts[0]->pos(), fedge->edge->verts[1]->pos(),
                v0, v1) < C_EDGES_INTERS_DIST * m_prefLen)
            return true;

    return false;
}


bool surface::Face::doesSegmentIntersectsWithFront(const pmg::Vert* v0, const vec3& v1) const
{
    for (auto& fedge : m_frontEdges)
    {
        if (fedge->edge->contains(v0))
            continue;

        if (spatalgs::segmentsDistance(
                fedge->edge->verts[0]->pos(), fedge->edge->verts[1]->pos(),
                v0->pos(), v1) < C_EDGES_INTERS_DIST * m_prefLen)
            return true;
    }

    return false;
}


vec3 surface::Face::computeNormalInTriang(sfront::Edge* fEdge, const vec3& oppVertPos)
{
    vec3 p0 = fEdge->edge->verts[0]->pos();
    vec3 p1 = fEdge->edge->verts[1]->pos();
    return (spatalgs::project(oppVertPos, p0, p1) - oppVertPos).normalize();
}


bool surface::Face::tryComputeNewVertPosType2(sfront::Edge* fEdge, vec3& out_pos)
{
    std::array<sfront::Vert*, 2> main_f_verts;
    main_f_verts[0] = findFrontVert(fEdge->edge->verts[0]);
    main_f_verts[1] = findFrontVert(fEdge->edge->verts[1]);

    auto adj_f_edges0 = main_f_verts[0]->findAdjEdges();
    auto adj_f_edges1 = main_f_verts[1]->findAdjEdges();
    auto en0 = std::get<0>(adj_f_edges0) == fEdge ? std::get<1>(adj_f_edges0) : std::get<0>(adj_f_edges0);
    auto en1 = std::get<0>(adj_f_edges1) == fEdge ? std::get<1>(adj_f_edges1) : std::get<0>(adj_f_edges1);

    vec3 vn0_pos = en0->edge->findNot(fEdge->edge)->pos();
    vec3 vn1_pos = en1->edge->findNot(fEdge->edge)->pos();

    auto v0 = main_f_verts[0]->vert;
    auto v1 = main_f_verts[1]->vert;

    vec3 e0 = (vn0_pos - 2.0 * v0->pos() + main_f_verts[1]->vert->pos()).normalize();
    vec3 e1 = (vn1_pos - 2.0 * v1->pos() + v0->pos()).normalize();

    vec3 new_pos = spatalgs::linesClosestPoint(v0->pos(), v0->pos() + e0, v1->pos(), v1->pos() + e1);

    vec3 v0_to_np = new_pos - v0->pos();
    vec3 v1_to_np = new_pos - v1->pos();
    if (doesSegmentIntersectsWithFront(v0, new_pos + K_MIN_DIS * v0_to_np) ||
        doesSegmentIntersectsWithFront(v1, new_pos + K_MIN_DIS * v1_to_np))
        return false;

    out_pos = new_pos;
    return true;
}


bool surface::Face::tryComputeNewVertPosType1(sfront::Edge* fEdge, vec3& out_pos, std::size_t smallAngleIdx)
{
    auto main_vert = fEdge->edge->verts[smallAngleIdx];
    auto sec_vert  = fEdge->edge->findNot(main_vert);
    auto main_f_vert = findFrontVert(main_vert);

    auto adj_f_edges = main_f_vert->findAdjEdges();
    auto en = std::get<0>(adj_f_edges) == fEdge ? std::get<1>(adj_f_edges) : std::get<0>(adj_f_edges);
    auto vn = en->edge->findNot(main_vert);

    vec3 e = (sec_vert->pos() - static_cast<real_t>(2.0) * main_vert->pos() + vn->pos()).normalize();

    real_t av_magn = static_cast<real_t>(0.5) * (fEdge->edge->magnitude() + en->edge->magnitude());
    real_t raw_deform = K_D * (m_prefLen - av_magn);
    real_t deform = raw_deform < av_magn * K_MAXD ? raw_deform : av_magn * K_MAXD;
    real_t magn_d = av_magn + deform;
    vec3 new_pos = main_vert->pos() + magn_d * e;

//    vec3 new_pos = main_vert->pos() + av_magn * e;

    vec3 v0_to_np = new_pos - main_vert->pos();
    vec3 v1_to_np = new_pos - sec_vert->pos();
    if (doesSegmentIntersectsWithFront(main_vert->pos(), new_pos + K_MIN_DIS * v0_to_np) ||
        doesSegmentIntersectsWithFront(sec_vert->pos(),  new_pos + K_MIN_DIS * v1_to_np))
        return false;

    out_pos = new_pos;
    return true;
}


bool surface::Face::tryComputeNewVertPosType0(sfront::Edge* fEdge, vec3& out_pos)
{
    real_t magn = fEdge->edge->magnitude();
    real_t raw_deform = K_D * (m_prefLen - magn);
    real_t deform = raw_deform < magn * K_MAXD ? raw_deform : magn * K_MAXD;
    real_t magn_d = magn + deform;
    vec3 new_pos = fEdge->computeCenter() + static_cast<real_t>(0.5) * std::sqrt(static_cast<real_t>(4.0) * magn_d * magn_d - magn * magn) * fEdge->normal;

//    vec3 new_pos = fEdge->computeCenter() + fEdge->normal * m_prefLen * SQRT3_2;

    auto v0 = fEdge->edge->verts[0];
    auto v1 = fEdge->edge->verts[1];
    vec3 v_to_np0 = (new_pos - v0->pos());
    vec3 v_to_np1 = (new_pos - v1->pos());
    if (doesSegmentIntersectsWithFront(v0, new_pos + v_to_np0 * K_MIN_DIS) ||
        doesSegmentIntersectsWithFront(v1, new_pos + v_to_np1 * K_MIN_DIS))
        return false;

    out_pos = new_pos;
    return true;
}


bool surface::Face::tryComputeNewVertPos(sfront::Edge* fEdge, vec3& out_pos)
{
    std::array<real_t, 2> angs
    {
        findFrontVert(fEdge->edge->verts[0])->angle(),
        findFrontVert(fEdge->edge->verts[1])->angle()
    };
    std::array<std::size_t, 2> idces;
    std::size_t n_small_angs = 0;
    if (angs[0] < degToRad(120)) idces[n_small_angs++] = 0;
    if (angs[1] < degToRad(120)) idces[n_small_angs++] = 1;

    switch (n_small_angs)
    {
    case 0: return tryComputeNewVertPosType0(fEdge, out_pos);
    case 1: return tryComputeNewVertPosType1(fEdge, out_pos, idces[0]);
    case 2: return tryComputeNewVertPosType2(fEdge, out_pos);
    }

    return false;
}


pair_rr surface::Face::computeMinMaxEdgesLengths(const vec3& p0, const vec3& p1, const vec3& p2)
{
    auto min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2);
    min_max.first  = std::sqrt(min_max.first);
    min_max.second = std::sqrt(min_max.second);
    return min_max;
}


pair_rr surface::Face::computeMinMaxEdgesSqrLengths(const vec3& p0, const vec3& p1, const vec3& p2)
{
    std::array<real_t, 3> sqr_magns;
    sqr_magns[0] = (p1 - p0).sqrMagnitude();
    sqr_magns[1] = (p2 - p0).sqrMagnitude();
    sqr_magns[2] = (p2 - p1).sqrMagnitude();
    return std::minmax({ sqr_magns[0], sqr_magns[1], sqr_magns[2] });
}


real_t surface::Face::computeTriangSimpleQuality(const vec3& p0, const vec3& p1, const vec3& p2)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2);
    return std::sqrt(sqr_min_max.first / sqr_min_max.second);
}


real_t surface::Face::computeTriangSimpleSqrQuality(const vec3& p0, const vec3& p1, const vec3& p2)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2);
    return sqr_min_max.first / sqr_min_max.second;
}


sfront::Edge* surface::Face::chooseEdgeForExhaustionWithNewVert(sfront::Vert* fVert)
{
    auto adj_edges = fVert->findAdjEdges();
    return std::get<0>(adj_edges)->edge->sqrMagnitude() < std::get<1>(adj_edges)->edge->sqrMagnitude() ?
                std::get<0>(adj_edges) : std::get<1>(adj_edges);
}


void surface::Face::exhaustWithNewVert(sfront::Edge* fEdge, const vec3& vertPos)
{
    std::array<pmg::Vert*, 2> main_verts = { fEdge->edge->verts[0], fEdge->edge->verts[1] };

    pmg::Vert* new_vert = addToFront(new pmg::Vert(vertPos))->vert;
    auto new_f_edge0 = addToFront(new pmg::Edge(main_verts[0], new_vert));
    auto new_f_edge1 = addToFront(new pmg::Edge(main_verts[1], new_vert));
    new_f_edge0->normal = computeNormalInTriang(new_f_edge0, fEdge->edge->findNot(new_f_edge0->edge)->pos());
    new_f_edge1->normal = computeNormalInTriang(new_f_edge1, fEdge->edge->findNot(new_f_edge1->edge)->pos());

    m_innerFaces.push_back(new pmg::Face(fEdge->edge, new_f_edge0->edge, new_f_edge1->edge));

    findFrontVert(main_verts[0])->refreshAngleData();
    findFrontVert(main_verts[1])->refreshAngleData();

    removeFromFront(fEdge);
}


void surface::Face::exhaustWithoutNewVert(sfront::Vert* fVert)
{
    auto adj_edges = fVert->findAdjEdges();
    std::pair<pmg::Vert*, pmg::Vert*> opp_verts =
    {
        adj_edges.first->edge->findNot(fVert->vert),
        adj_edges.second->edge->findNot(fVert->vert)
    };

    auto new_f_edge = addToFront(new pmg::Edge(opp_verts.first, opp_verts.second));
    new_f_edge->normal = computeNormalInTriang(new_f_edge, fVert->vert->pos());

    m_innerFaces.push_back(new pmg::Face(adj_edges.first->edge, adj_edges.second->edge, new_f_edge->edge));

    findFrontVert(opp_verts.first)->refreshAngleData();
    findFrontVert(opp_verts.second)->refreshAngleData();

    removeFromFront(fVert);
    removeFromFront(adj_edges.first);
    removeFromFront(adj_edges.second);
}


bool surface::Face::tryExhaustWithoutNewVert(sfront::Vert* fVert)
{
    if (anyVertInsidePotentialTriangCheck(fVert))
        return false;

    exhaustWithoutNewVert(fVert);
    return true;
}


bool surface::Face::tryExhaustWithNewVert(sfront::Vert* fVert)
{
    auto exhaust_f_edge = chooseEdgeForExhaustionWithNewVert(fVert);
    vec3 new_vert_pos;
    if (!tryComputeNewVertPos(exhaust_f_edge, new_vert_pos))
        return false;

    exhaustWithNewVert(exhaust_f_edge, new_vert_pos);
    return true;
}


sfront::Vert* surface::Face::currentFrontVert(real_t maxCompl) const
{
    real_t cur_max_compl = 0.0;
    sfront::Vert* cur_max_f_edge = nullptr;
    for (auto& f_vert : m_frontVerts)
    {
        real_t cur_compl = f_vert->complexity();
        if (   cur_compl > cur_max_compl
            && cur_compl < maxCompl)
        {
            cur_max_compl  = cur_compl;
            cur_max_f_edge = f_vert;
        }
    }

    return cur_max_f_edge;
}


bool surface::Face::exhaustWithoutNewVertPriorityPredicate(sfront::Vert* fEdge)
{
    if (fEdge->angle() < degToRad(60))
        return true;

    if (fEdge->angle() > degToRad(80))
        return false;

    auto adj_edges = fEdge->findAdjEdges();
    real_t div = adj_edges.first->edge->sqrMagnitude() / adj_edges.second->edge->sqrMagnitude();
    if (NINE_DIV_SIXTEEN < div && div < SIXTEEN_DIV_NINE)
        return true;

    return false;
}


bool surface::Face::exhaustWithNewVertPriorityPredicate(sfront::Vert* fEdge)
{
    if (fEdge->angle() > degToRad(100))
        return true;

    return false;
}


surface::Face::ExhaustType surface::Face::computeExhaustionTypeQualityPriority(
        sfront::Vert* fVert, sfront::Edge*& out_withNWFrontEdge, vec3*& out_withNWNewVertPos)
{
    if (anyVertInsidePotentialTriangCheck(fVert))
        return ExhaustType::WithNewVert;

    auto opp_verts = fVert->oppVerts();
    real_t without_nv_quality = computeTriangSimpleSqrQuality(
        fVert->vert->pos(),
        std::get<0>(opp_verts)->pos(),
        std::get<1>(opp_verts)->pos());

    sfront::Edge* f_edge = chooseEdgeForExhaustionWithNewVert(fVert);
    vec3 new_vert_pos;
    if (!tryComputeNewVertPos(f_edge, new_vert_pos))
        return ExhaustType::DontExhaust;

    real_t with_nv_quality = computeTriangSimpleSqrQuality(
        f_edge->edge->verts[0]->pos(),
        f_edge->edge->verts[1]->pos(),
        new_vert_pos);

    if (without_nv_quality > with_nv_quality)
        return ExhaustType::WithoutNewVert;

    out_withNWFrontEdge = f_edge;
    out_withNWNewVertPos = new vec3(new_vert_pos);
    return ExhaustType::WithNewVert;
}


void surface::Face::processLastFace()
{
    std::array<pmg::Edge*, 3> edges;
    std::size_t i = 0;
    for (auto& fedge : m_frontEdges)
        edges[i++] = fedge->edge;

    m_innerFaces.push_back(new pmg::Face(edges[0], edges[1], edges[2]));

    for (auto& fedge : m_frontEdges)
        delete fedge;
    m_frontEdges.clear();

    for (auto& fvert : m_frontVerts)
        delete fvert;
    m_frontVerts.clear();
}


void surface::Face::processAngles()
{
    if (m_frontEdges.size() < 3)
        throw std::logic_error("Wrong input data.\nError in function: pmg::surface::Face::processAngles");

    if (m_frontEdges.size() == 3)
    {
        processLastFace();
        return;
    }

    real_t max_compl = std::numeric_limits<real_t>::max();
    for (sfront::Vert* cur_f_vert = currentFrontVert(max_compl);; cur_f_vert = currentFrontVert(max_compl))
    {
        if (!cur_f_vert)
            throw std::logic_error("pmg::surface::Face::currentFrontVert returned nullptr");

        if (exhaustWithoutNewVertPriorityPredicate(cur_f_vert))
        {
            if (!tryExhaustWithoutNewVert(cur_f_vert))
            {
                max_compl = cur_f_vert->complexity();
                continue;
            }
        }
        else if (exhaustWithNewVertPriorityPredicate(cur_f_vert))
        {
            if (!tryExhaustWithNewVert(cur_f_vert))
            {
                max_compl = cur_f_vert->complexity();
                continue;
            }
        }
        else
        {
            sfront::Edge* exhaust_from_f_edge = nullptr;
            vec3* new_vert_pos = nullptr;
            switch (computeExhaustionTypeQualityPriority(cur_f_vert, exhaust_from_f_edge, new_vert_pos))
            {
            case ExhaustType::WithoutNewVert:
                exhaustWithoutNewVert(cur_f_vert);
                break;

            case ExhaustType::WithNewVert:
                if (new_vert_pos)
                {
                    exhaustWithNewVert(exhaust_from_f_edge, *new_vert_pos);
                    delete new_vert_pos;
                }
                else
                {
                    if (!tryExhaustWithNewVert(cur_f_vert))
                    {
                        max_compl = cur_f_vert->complexity();
                        continue;
                    }
                }
                break;

            case ExhaustType::DontExhaust:
                max_compl = cur_f_vert->complexity();
                continue;
            }
        }
        max_compl = std::numeric_limits<real_t>::max();

        if (m_frontEdges.size() == 3)
        {
            processLastFace();
            return;
        }
    }
}


void surface::Face::smoothMesh(std::size_t nIters)
{
    for (std::size_t i = 0; i < nIters; i++)
    {
        for (auto& vert : m_innerVerts)
        {
            vec3 shift;
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


pair_ff surface::Face::find2AdjFaces(pmg::Edge* edge) const
{
    pair_ff res;
    bool not_found_yet = true;
    for (auto& face : m_innerFaces)
    {
        if (face->contains(edge))
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

    throw std::logic_error("pmg::surface::Face::find2AdjFaces didn't find 2 adjacent Faces.");
}


bool surface::Face::flipIfNeeded(pmg::Edge* edge)
{
    std::array<pmg::Vert*, 2> opp_nodes;
    auto around_faces = find2AdjFaces(edge);
    opp_nodes[0] = std::get<0>(around_faces)->findVertNot(edge);
    opp_nodes[1] = std::get<1>(around_faces)->findVertNot(edge);

    real_t alpha = std::acos(vec3::cos(*edge->verts[0] - *opp_nodes[0], *edge->verts[1] - *opp_nodes[0]));
    real_t beta  = std::acos(vec3::cos(*edge->verts[0] - *opp_nodes[1], *edge->verts[1] - *opp_nodes[1]));

    if (alpha + beta <= PI)
        return false;


    auto old_edge   = std::find(m_innerEdges.begin(), m_innerEdges.end(), edge);
    auto old_face0 = std::find(m_innerFaces.begin(), m_innerFaces.end(), std::get<0>(around_faces));
    auto old_face1 = std::find(m_innerFaces.begin(), m_innerFaces.end(), std::get<1>(around_faces));

    pmg::Edge* new_edge = new pmg::Edge(opp_nodes[0], opp_nodes[1]);
    pmg::Face* new_face0 = new pmg::Face(
            std::get<0>(around_faces)->findEdge(opp_nodes[0], edge->verts[0]),
            std::get<1>(around_faces)->findEdge(opp_nodes[1], edge->verts[0]),
            new_edge);
    pmg::Face* new_face1 = new pmg::Face(
            std::get<0>(around_faces)->findEdge(opp_nodes[0], edge->verts[1]),
            std::get<1>(around_faces)->findEdge(opp_nodes[1], edge->verts[1]),
            new_edge);

    delete *old_edge;
    delete *old_face0;
    delete *old_face1;

    *old_edge = new_edge;
    *old_face0 = new_face0;
    *old_face1 = new_face1;

    return true;
}


void surface::Face::delaunayPostP()
{
    for (auto& edge : m_innerEdges)
        flipIfNeeded(edge);
}


void surface::Face::optimizeMesh(std::size_t nSmoothIters, std::size_t nDelaunaySmoothIters)
{
    for (std::size_t i = 0; i < nDelaunaySmoothIters; i++)
    {
        delaunayPostP();
        smoothMesh(nSmoothIters);
    }
}


void surface::Face::computeFrontNormals() const
{
    for (auto& fedge : m_frontEdges)
        fedge->computeNormal();
}


void surface::Face::initializeFront()
{
    std::vector<surface::Vert*> sverts;
    for (auto& this_edge : edges)
        for (auto& svert : this_edge->verts)
            if (std::find(sverts.begin(), sverts.end(), svert) == sverts.end())
                sverts.push_back(svert);

    for (auto& svert : sverts)
        m_frontVerts.push_back(new sfront::Vert(this, svert->attachedVert));
    sverts.clear();

    for (auto& this_edge : edges)
        for (auto& vert : this_edge->innerVerts())
            m_frontVerts.push_back(new sfront::Vert(this, vert));

    for (auto& this_edge : edges)
        for (auto& edge : this_edge->innerEdges())
            m_frontEdges.push_back(new sfront::Edge(this, edge));
}
