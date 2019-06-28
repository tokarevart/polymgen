// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "face.h"
#include <algorithm>
#include <cmath>
#include "../../helpers/mathconsts.h"
#include "../../helpers/spatial/algs.h"


using namespace pmg;
namespace sfront = surface::front;
using pair_rr = std::pair<real_t, real_t>;
using pair_ff = std::pair<pmg::Face*, pmg::Face*>;
using pair_ee = std::pair<pmg::Edge*, pmg::Edge*>;
using vec3 = spt::vec<3, real_t>;


#define NINE_DIV_SIXTEEN static_cast<real_t>(0.5625)
#define SIXTEEN_DIV_NINE static_cast<real_t>(1.7777777777777777)
#define SQRT3_2          static_cast<real_t>(0.8660254037844386)

#define K_MIN_DIS              static_cast<real_t>(2e-1)
#define C_EDGES_INTERS_DIST static_cast<real_t>(1e-4)

#define K_MAXD static_cast<real_t>(0.3)
#define K_D    static_cast<real_t>(0.5)


template <typename T>
constexpr real_t degToRad(T value) {
    return value * DEG_IN_RAD;
}


real_t surface::Face::preferred_length() const {
    return m_prefLen;
}


const std::list<pmg::Face*>& surface::Face::inner_faces() const {
    return m_inner_faces;
}


const std::list<pmg::Edge*>& surface::Face::inner_edges() const {
    return m_inner_edges;
}


const std::list<pmg::Vert*>& surface::Face::inner_verts() const {
    return m_inner_verts;
}


const std::list<sfront::Edge*>& surface::Face::front_edges() const {
    return m_front_edges;
}


surface::Vert* surface::Face::find_vert_not(const surface::Edge* edge) const {
    for (auto& face_edge : edges) {
        if (edge != face_edge) {
            if (!edge->contains(face_edge->verts[0]))
                return face_edge->verts[0];
            else
                return face_edge->verts[1];
        }
    }

    return nullptr;
}


surface::Edge* surface::Face::find_surface_edge_containing(const pmg::Edge* edge) const {
    for (auto& sedge : edges)
        if (sedge->contains(edge))
            return sedge;

    return nullptr;
}


void surface::Face::triangulate(real_t preferredLen, genparams::Surface gen_params) {
    m_prefLen = preferredLen;
    initializeFront();
    computeFrontNormals();
    processAngles();
//    if (globalIntersectionCheck())
//        throw std::logic_error("Intersection error.\npmg::surface::Face::globalIntersectionCheck returned true.");

    optimize_mesh(gen_params.nSmoothIters, gen_params.nDelaunaySmoothIters);
}


bool surface::Face::contains(const surface::Edge* edge) const {
    for (auto& edge_ : edges)
        if (edge_ == edge)
            return true;

    return false;
}


bool surface::Face::contains(const surface::Vert* vert) const {
    for (auto& edge : edges)
        if (edge->contains(vert))
            return true;

    return false;
}




surface::Face::Face(
    const surface::Edge* edge0,
    const surface::Edge* edge1,
    const surface::Edge* edge2) {
    edges[0] = const_cast<surface::Edge*>(edge0);
    edges[1] = const_cast<surface::Edge*>(edge1);
    edges[2] = const_cast<surface::Edge*>(edge2);
}




sfront::Vert* surface::Face::find_front_vert(const pmg::Vert* vert) const {
    for (auto& fvert : m_front_verts)
        if (fvert->vert == vert)
            return fvert;

    return nullptr;
}


sfront::Edge* surface::Face::add_to_front(const pmg::Edge* edge) {
    sfront::Edge* new_f_edge = new sfront::Edge(this, edge);
    m_front_edges.push_back(new_f_edge);
    m_inner_edges.push_back(new_f_edge->edge);
    return new_f_edge;
}


sfront::Vert* surface::Face::add_to_front(const pmg::Vert* vert) {
    sfront::Vert* new_f_vert = new sfront::Vert(this, vert);
    m_front_verts.push_back(new_f_vert);
    m_inner_verts.push_back(new_f_vert->vert);
    return new_f_vert;
}


void surface::Face::remove_from_front(sfront::Edge* front_edge) {
    m_front_edges.erase(std::find(m_front_edges.begin(), m_front_edges.end(), front_edge));
    delete front_edge;
}


void surface::Face::remove_from_front(sfront::Vert* fVert) {
    m_front_verts.erase(std::find(m_front_verts.begin(), m_front_verts.end(), fVert));
    delete fVert;
}


bool surface::Face::any_vert_inside_potential_triangle_check(sfront::Vert* fVert) const {
    auto opp_verts_l = fVert->opp_verts();
    std::array<pmg::Vert*, 3> tr
    {
        fVert->vert,
        opp_verts_l.first,
        opp_verts_l.second
    };

    for (auto& fvert : m_front_verts)
        if (fvert->vert != tr[0] && fvert->vert != tr[1] && fvert->vert != tr[2] &&
            spt::is_point_on_triangle(fvert->vert->pos(), tr[0]->pos(), tr[1]->pos(), tr[2]->pos()))
            return true;

    return false;
}


bool surface::Face::does_segment_intersects_with_front(const vec3& v0, const vec3& v1) const {
    for (auto& f_edge : m_front_edges)
        if (spt::segments_distance(
            f_edge->edge->verts[0]->pos(), f_edge->edge->verts[1]->pos(),
            v0, v1) < C_EDGES_INTERS_DIST * m_prefLen)
            return true;

    return false;
}


bool surface::Face::does_segment_intersects_with_front(const pmg::Vert* v0, const vec3& v1) const {
    for (auto& f_edge : m_front_edges) {
        if (f_edge->edge->contains(v0))
            continue;

        if (spt::segments_distance(
            f_edge->edge->verts[0]->pos(), f_edge->edge->verts[1]->pos(),
            v0->pos(), v1) < C_EDGES_INTERS_DIST * m_prefLen)
            return true;
    }

    return false;
}


vec3 surface::Face::normal_in_triangle(sfront::Edge* front_edge, const vec3& opp_vert_pos) {
    vec3 p0 = front_edge->edge->verts[0]->pos();
    vec3 p1 = front_edge->edge->verts[1]->pos();
    return (spt::project(opp_vert_pos, p0, p1) - opp_vert_pos).normalize();
}


bool surface::Face::tryComputeNewVertPosType2(sfront::Edge* front_edge, vec3& out_pos) {
    std::array<sfront::Vert*, 2> main_f_verts;
    main_f_verts[0] = find_front_vert(front_edge->edge->verts[0]);
    main_f_verts[1] = find_front_vert(front_edge->edge->verts[1]);

    auto adj_f_edges0 = main_f_verts[0]->findAdjEdges();
    auto adj_f_edges1 = main_f_verts[1]->findAdjEdges();
    auto en0 = std::get<0>(adj_f_edges0) == front_edge ? std::get<1>(adj_f_edges0) : std::get<0>(adj_f_edges0);
    auto en1 = std::get<0>(adj_f_edges1) == front_edge ? std::get<1>(adj_f_edges1) : std::get<0>(adj_f_edges1);

    vec3 vn0_pos = en0->edge->findNot(front_edge->edge)->pos();
    vec3 vn1_pos = en1->edge->findNot(front_edge->edge)->pos();

    auto v0 = main_f_verts[0]->vert;
    auto v1 = main_f_verts[1]->vert;

    vec3 e0 = (vn0_pos - v0->pos() * 2.0 + main_f_verts[1]->vert->pos()).normalize();
    vec3 e1 = (vn1_pos - v1->pos() * 2.0 + v0->pos()).normalize();

    vec3 new_pos = spt::lines_closest_point(v0->pos(), v0->pos() + e0, v1->pos(), v1->pos() + e1);

    vec3 v0_to_np = new_pos - v0->pos();
    vec3 v1_to_np = new_pos - v1->pos();
    if (does_segment_intersects_with_front(v0, new_pos + v0_to_np * K_MIN_DIS) ||
        does_segment_intersects_with_front(v1, new_pos + v1_to_np * K_MIN_DIS))
        return false;

    out_pos = new_pos;
    return true;
}


bool surface::Face::tryComputeNewVertPosType1(sfront::Edge* front_edge, vec3& out_pos, std::size_t smallAngleIdx) {
    auto main_vert = front_edge->edge->verts[smallAngleIdx];
    auto sec_vert = front_edge->edge->findNot(main_vert);
    auto main_f_vert = find_front_vert(main_vert);

    auto adj_f_edges = main_f_vert->findAdjEdges();
    auto en = std::get<0>(adj_f_edges) == front_edge ? std::get<1>(adj_f_edges) : std::get<0>(adj_f_edges);
    auto vn = en->edge->findNot(main_vert);

    vec3 e = (sec_vert->pos() - main_vert->pos() * static_cast<real_t>(2.0) + vn->pos()).normalize();

    real_t av_magn = static_cast<real_t>(0.5) * (front_edge->edge->magnitude() + en->edge->magnitude());
    real_t raw_deform = K_D * (m_prefLen - av_magn);
    real_t deform = raw_deform < av_magn * K_MAXD ? raw_deform : av_magn * K_MAXD;
    real_t magn_d = av_magn + deform;
    vec3 new_pos = main_vert->pos() + e * magn_d;

//    vec3 new_pos = main_vert->pos() + av_magn * e;

    vec3 v0_to_np = new_pos - main_vert->pos();
    vec3 v1_to_np = new_pos - sec_vert->pos();
    if (does_segment_intersects_with_front(main_vert->pos(), new_pos + v0_to_np * K_MIN_DIS) ||
        does_segment_intersects_with_front(sec_vert->pos(), new_pos + v1_to_np * K_MIN_DIS))
        return false;

    out_pos = new_pos;
    return true;
}


bool surface::Face::tryComputeNewVertPosType0(sfront::Edge* front_edge, vec3& out_pos) {
    real_t magn = front_edge->edge->magnitude();
    real_t raw_deform = K_D * (m_prefLen - magn);
    real_t deform = raw_deform < magn * K_MAXD ? raw_deform : magn * K_MAXD;
    real_t magn_d = magn + deform;
    vec3 new_pos = front_edge->center() + front_edge->normal * (static_cast<real_t>(0.5) * std::sqrt(static_cast<real_t>(4.0) * magn_d * magn_d - magn * magn));

//    vec3 new_pos = front_edge->center() + front_edge->normal * m_prefLen * SQRT3_2;

    auto v0 = front_edge->edge->verts[0];
    auto v1 = front_edge->edge->verts[1];
    vec3 v_to_np0 = (new_pos - v0->pos());
    vec3 v_to_np1 = (new_pos - v1->pos());
    if (does_segment_intersects_with_front(v0, new_pos + v_to_np0 * K_MIN_DIS) ||
        does_segment_intersects_with_front(v1, new_pos + v_to_np1 * K_MIN_DIS))
        return false;

    out_pos = new_pos;
    return true;
}


bool surface::Face::tryComputeNewVertPos(sfront::Edge* front_edge, vec3& out_pos) {
    std::array<real_t, 2> angs
    {
        find_front_vert(front_edge->edge->verts[0])->angle(),
        find_front_vert(front_edge->edge->verts[1])->angle()
    };
    std::array<std::size_t, 2> idces;
    std::size_t n_small_angs = 0;
    if (angs[0] < degToRad(120)) idces[n_small_angs++] = 0;
    if (angs[1] < degToRad(120)) idces[n_small_angs++] = 1;

    switch (n_small_angs) {
    case 0: return tryComputeNewVertPosType0(front_edge, out_pos);
    case 1: return tryComputeNewVertPosType1(front_edge, out_pos, idces[0]);
    case 2: return tryComputeNewVertPosType2(front_edge, out_pos);
    }

    return false;
}


pair_rr surface::Face::computeMinMaxEdgesLengths(const vec3& p0, const vec3& p1, const vec3& p2) {
    auto min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2);
    min_max.first = std::sqrt(min_max.first);
    min_max.second = std::sqrt(min_max.second);
    return min_max;
}


pair_rr surface::Face::computeMinMaxEdgesSqrLengths(const vec3& p0, const vec3& p1, const vec3& p2) {
    std::array<real_t, 3> sqr_magns;
    sqr_magns[0] = (p1 - p0).sqr_magnitude();
    sqr_magns[1] = (p2 - p0).sqr_magnitude();
    sqr_magns[2] = (p2 - p1).sqr_magnitude();
    return std::minmax({ sqr_magns[0], sqr_magns[1], sqr_magns[2] });
}


real_t surface::Face::computeTriangSimpleQuality(const vec3& p0, const vec3& p1, const vec3& p2) {
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2);
    return std::sqrt(sqr_min_max.first / sqr_min_max.second);
}


real_t surface::Face::computeTriangSimpleSqrQuality(const vec3& p0, const vec3& p1, const vec3& p2) {
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2);
    return sqr_min_max.first / sqr_min_max.second;
}


sfront::Edge* surface::Face::chooseEdgeForExhaustionWithNewVert(sfront::Vert* fVert) {
    auto adj_edges = fVert->findAdjEdges();
    return std::get<0>(adj_edges)->edge->sqr_magnitude() < std::get<1>(adj_edges)->edge->sqr_magnitude() ?
        std::get<0>(adj_edges) : std::get<1>(adj_edges);
}


void surface::Face::exhaustWithNewVert(sfront::Edge* front_edge, const vec3& vertPos) {
    std::array<pmg::Vert*, 2> main_verts = { front_edge->edge->verts[0], front_edge->edge->verts[1] };

    pmg::Vert* new_vert = add_to_front(new pmg::Vert(vertPos))->vert;
    auto new_f_edge0 = add_to_front(new pmg::Edge(main_verts[0], new_vert));
    auto new_f_edge1 = add_to_front(new pmg::Edge(main_verts[1], new_vert));
    new_f_edge0->normal = normal_in_triangle(new_f_edge0, front_edge->edge->findNot(new_f_edge0->edge)->pos());
    new_f_edge1->normal = normal_in_triangle(new_f_edge1, front_edge->edge->findNot(new_f_edge1->edge)->pos());

    m_inner_faces.push_back(new pmg::Face(front_edge->edge, new_f_edge0->edge, new_f_edge1->edge));

    find_front_vert(main_verts[0])->refresh_angle_data();
    find_front_vert(main_verts[1])->refresh_angle_data();

    remove_from_front(front_edge);
}


void surface::Face::exhaustWithoutNewVert(sfront::Vert* fVert) {
    auto adj_edges = fVert->findAdjEdges();
    std::pair<pmg::Vert*, pmg::Vert*> opp_verts_l =
    {
        adj_edges.first->edge->findNot(fVert->vert),
        adj_edges.second->edge->findNot(fVert->vert)
    };

    auto new_f_edge = add_to_front(new pmg::Edge(opp_verts_l.first, opp_verts_l.second));
    new_f_edge->normal = normal_in_triangle(new_f_edge, fVert->vert->pos());

    m_inner_faces.push_back(new pmg::Face(adj_edges.first->edge, adj_edges.second->edge, new_f_edge->edge));

    find_front_vert(opp_verts_l.first)->refresh_angle_data();
    find_front_vert(opp_verts_l.second)->refresh_angle_data();

    remove_from_front(fVert);
    remove_from_front(adj_edges.first);
    remove_from_front(adj_edges.second);
}


bool surface::Face::tryExhaustWithoutNewVert(sfront::Vert* fVert) {
    if (any_vert_inside_potential_triangle_check(fVert))
        return false;

    exhaustWithoutNewVert(fVert);
    return true;
}


bool surface::Face::tryExhaustWithNewVert(sfront::Vert* fVert) {
    auto exhaust_f_edge = chooseEdgeForExhaustionWithNewVert(fVert);
    vec3 new_vert_pos;
    if (!tryComputeNewVertPos(exhaust_f_edge, new_vert_pos))
        return false;

    exhaustWithNewVert(exhaust_f_edge, new_vert_pos);
    return true;
}


sfront::Vert* surface::Face::currentFrontVert(real_t maxCompl) const {
    real_t cur_max_compl = 0.0;
    sfront::Vert* cur_max_f_edge = nullptr;
    for (auto& f_vert : m_front_verts) {
        real_t cur_compl = f_vert->complexity();
        if (cur_compl > cur_max_compl
            && cur_compl < maxCompl) {
            cur_max_compl = cur_compl;
            cur_max_f_edge = f_vert;
        }
    }

    return cur_max_f_edge;
}


bool surface::Face::exhaustWithoutNewVertPriorityPredicate(sfront::Vert* front_edge) {
    if (front_edge->angle() < degToRad(60))
        return true;

    if (front_edge->angle() > degToRad(80))
        return false;

    auto adj_edges = front_edge->findAdjEdges();
    real_t div = adj_edges.first->edge->sqr_magnitude() / adj_edges.second->edge->sqr_magnitude();
    if (NINE_DIV_SIXTEEN < div && div < SIXTEEN_DIV_NINE)
        return true;

    return false;
}


bool surface::Face::exhaustWithNewVertPriorityPredicate(sfront::Vert* front_edge) {
    if (front_edge->angle() > degToRad(100))
        return true;

    return false;
}


surface::Face::ExhaustType surface::Face::computeExhaustionTypeQualityPriority(
    sfront::Vert* fVert, sfront::Edge*& out_withNWFrontEdge, vec3*& out_withNWNewVertPos) {
    if (any_vert_inside_potential_triangle_check(fVert))
        return ExhaustType::WithNewVert;

    auto opp_verts_l = fVert->opp_verts();
    real_t without_nv_quality = computeTriangSimpleSqrQuality(
        fVert->vert->pos(),
        std::get<0>(opp_verts_l)->pos(),
        std::get<1>(opp_verts_l)->pos());

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


void surface::Face::processLastFace() {
    std::array<pmg::Edge*, 3> edges;
    std::size_t i = 0;
    for (auto& f_edge : m_front_edges)
        edges[i++] = f_edge->edge;

    m_inner_faces.push_back(new pmg::Face(edges[0], edges[1], edges[2]));

    for (auto& f_edge : m_front_edges)
        delete f_edge;
    m_front_edges.clear();

    for (auto& fvert : m_front_verts)
        delete fvert;
    m_front_verts.clear();
}


void surface::Face::processAngles() {
    if (m_front_edges.size() < 3)
        throw std::logic_error("Wrong input data.\nError in function: pmg::surface::Face::processAngles");

    if (m_front_edges.size() == 3) {
        processLastFace();
        return;
    }

    real_t max_compl = std::numeric_limits<real_t>::max();
    for (sfront::Vert* cur_f_vert = currentFrontVert(max_compl);; cur_f_vert = currentFrontVert(max_compl)) {
        if (!cur_f_vert)
            throw std::logic_error("pmg::surface::Face::currentFrontVert returned nullptr");

        if (exhaustWithoutNewVertPriorityPredicate(cur_f_vert)) {
            if (!tryExhaustWithoutNewVert(cur_f_vert)) {
                max_compl = cur_f_vert->complexity();
                continue;
            }
        } else if (exhaustWithNewVertPriorityPredicate(cur_f_vert)) {
            if (!tryExhaustWithNewVert(cur_f_vert)) {
                max_compl = cur_f_vert->complexity();
                continue;
            }
        } else {
            sfront::Edge* exhaust_from_f_edge = nullptr;
            vec3* new_vert_pos = nullptr;
            switch (computeExhaustionTypeQualityPriority(cur_f_vert, exhaust_from_f_edge, new_vert_pos)) {
            case ExhaustType::WithoutNewVert:
                exhaustWithoutNewVert(cur_f_vert);
                break;

            case ExhaustType::WithNewVert:
                if (new_vert_pos) {
                    exhaustWithNewVert(exhaust_from_f_edge, *new_vert_pos);
                    delete new_vert_pos;
                } else {
                    if (!tryExhaustWithNewVert(cur_f_vert)) {
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

        if (m_front_edges.size() == 3) {
            processLastFace();
            return;
        }
    }
}


void surface::Face::smooth_mesh(std::size_t nIters) {
    for (std::size_t i = 0; i < nIters; i++) {
        for (auto& vert : m_inner_verts) {
            vec3 shift;
            int delta_shifts_num = 0;
            for (auto& edge : m_inner_edges) {
                if (vert == edge->verts[0]) {
                    shift += edge->verts[1]->pos() - vert->pos();
                    delta_shifts_num++;
                } else if (vert == edge->verts[1]) {
                    shift += edge->verts[0]->pos() - vert->pos();
                    delta_shifts_num++;
                }
            }
            shift /= delta_shifts_num;
            vert->pos() = vert->pos() + shift;
        }
    }
}


pair_ff surface::Face::find2AdjFaces(pmg::Edge* edge) const {
    pair_ff res;
    bool not_found_yet = true;
    for (auto& face : m_inner_faces) {
        if (face->contains(edge)) {
            if (not_found_yet) {
                res.first = face;
                not_found_yet = false;
            } else {
                res.second = face;
                return res;
            }
        }
    }

    throw std::logic_error("pmg::surface::Face::find2AdjFaces didn't find 2 adjacent Faces.");
}


bool surface::Face::flipIfNeeded(pmg::Edge* edge) {
    std::array<pmg::Vert*, 2> opp_nodes;
    auto around_faces = find2AdjFaces(edge);
    opp_nodes[0] = std::get<0>(around_faces)->find_vert_not(edge);
    opp_nodes[1] = std::get<1>(around_faces)->find_vert_not(edge);

    real_t alpha = std::acos(spt::cos(edge->verts[0]->pos() - opp_nodes[0]->pos(), edge->verts[1]->pos() - opp_nodes[0]->pos()));
    real_t beta = std::acos(spt::cos(edge->verts[0]->pos() - opp_nodes[1]->pos(), edge->verts[1]->pos() - opp_nodes[1]->pos()));

    if (alpha + beta <= PI)
        return false;


    auto old_edge = std::find(m_inner_edges.begin(), m_inner_edges.end(), edge);
    auto old_face0 = std::find(m_inner_faces.begin(), m_inner_faces.end(), std::get<0>(around_faces));
    auto old_face1 = std::find(m_inner_faces.begin(), m_inner_faces.end(), std::get<1>(around_faces));

    pmg::Edge* new_edge = new pmg::Edge(opp_nodes[0], opp_nodes[1]);
    pmg::Face* new_face0 = new pmg::Face(
        std::get<0>(around_faces)->findEdge(opp_nodes[0], edge->verts[0]),
        std::get<1>(around_faces)->findEdge(opp_nodes[1], edge->verts[0]),
        new_edge);
    pmg::Face* new_face1 = new pmg::Face(
        std::get<0>(around_faces)->findEdge(opp_nodes[0], edge->verts[1]),
        std::get<1>(around_faces)->findEdge(opp_nodes[1], edge->verts[1]),
        new_edge);

    delete* old_edge;
    delete* old_face0;
    delete* old_face1;

    *old_edge = new_edge;
    *old_face0 = new_face0;
    *old_face1 = new_face1;

    return true;
}


void surface::Face::delaunay_postp() {
    for (auto& edge : m_inner_edges)
        flipIfNeeded(edge);
}


void surface::Face::optimize_mesh(std::size_t nSmoothIters, std::size_t nDelaunaySmoothIters) {
    for (std::size_t i = 0; i < nDelaunaySmoothIters; i++) {
        delaunay_postp();
        smooth_mesh(nSmoothIters);
    }
}


void surface::Face::computeFrontNormals() const {
    for (auto& f_edge : m_front_edges)
        f_edge->compute_normal();
}


void surface::Face::initializeFront() {
    std::vector<surface::Vert*> sverts;
    for (auto& this_edge : edges)
        for (auto& svert : this_edge->verts)
            if (std::find(sverts.begin(), sverts.end(), svert) == sverts.end())
                sverts.push_back(svert);

    for (auto& svert : sverts)
        m_front_verts.push_back(new sfront::Vert(this, svert->attached_vert));
    sverts.clear();

    for (auto& this_edge : edges)
        for (auto& vert : this_edge->inner_verts())
            m_front_verts.push_back(new sfront::Vert(this, vert));

    for (auto& this_edge : edges)
        for (auto& edge : this_edge->inner_edges())
            m_front_edges.push_back(new sfront::Edge(this, edge));
}
