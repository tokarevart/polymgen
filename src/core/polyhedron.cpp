// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "polyhedron.h"
#include <algorithm>
#include <iostream>
#include "../helpers/mathconsts.h"
#include "../helpers/spatial/algs.h"


using namespace pmg;
using pair_rr = std::pair<real_t, real_t>;
using pair_ff = std::pair<front::Face*, front::Face*>;
using vec3 = spt::vec<3, real_t>;


#define DET(a, b, c, d) \
        (a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0[0], corner1[0], point[0]) && \
         BETWEEN(corner0[1], corner1[1], point[1]))


// TODO: reduce the amount of defines
#define ALPHA_P      static_cast<real_t>(70.52877936550931)

#define ONE_3                static_cast<real_t>(0.3333333333333333)
#define SQRT3_2              static_cast<real_t>(0.8660254037844386)
#define SQRT_2_3             static_cast<real_t>(0.8164965809277260)
#define ONE_PLUS_SQRT2_SQRT3 static_cast<real_t>(1.3938468501173517)

#define C_MIN_DIS           static_cast<real_t>(2e-1)
#define C_EDGES_INTERS_DIST static_cast<real_t>(4e-3)

#define C_MAXD static_cast<real_t>(0.3)
#define C_D    static_cast<real_t>(0.4)


template <typename T>
constexpr real_t degToRad(T value) {
    return value * DEG_IN_RAD;
}


#define DEBUG


void Polyhedron::initialize_fface_fedges(front::Face* fface) const {
    int n_added = 0;
    for (auto& l_fedge : m_front_edges) {
        if (fface->x->contains(l_fedge->x)) {
            fface->add_front_edge(l_fedge);
            if (++n_added == 3)
                return;
        }
    }

    throw std::logic_error("Polyhedron::initialize_fface_fedges didn't find 3 front edges.");
}


void Polyhedron::initialize_front() {
    for (auto& svert : m_shell_verts)
        m_front_verts.push_back(new front::Vert(svert->attached_vert));

    for (auto& sedge : m_shell_edges)
        for (auto& vert : sedge->inner_verts())
            m_front_verts.push_back(new front::Vert(vert));

    for (auto& sface : m_shell_faces)
        for (auto& vert : sface->inner_verts())
            m_front_verts.push_back(new front::Vert(vert));

    for (auto& sedge : m_shell_edges)
        for (auto& edge : sedge->inner_edges())
            m_front_edges.push_back(new front::Edge(this, edge));

    for (auto& sface : m_shell_faces)
        for (auto& edge : sface->inner_edges())
            m_front_edges.push_back(new front::Edge(this, edge));

    for (auto& sface : m_shell_faces)
        for (auto& face : sface->inner_faces())
            m_front_faces.push_back(new front::Face(this, face));

    for (auto& fedge : m_front_edges)
        initialize_fedge_fverts(fedge);

    for (auto& fface : m_front_faces)
        initialize_fface_fedges(fface);
}


void Polyhedron::compute_front_normals() {
    for (auto& face : m_front_faces)
        face->compute_normal();
}


void pmg::Polyhedron::initialize_fedge_fverts(front::Edge* fedge) const {
    int n_added = 0;
    for (auto& l_fvert : m_front_verts)
        if (fedge->x->contains(l_fvert->x)) {
            fedge->add_front_vert(l_fvert);
            if (++n_added == 2)
                return;
        }

    throw std::logic_error("Polyhedron::initialize_fedge_fverts didn't find 2 front verts.");
}


shell::Edge* Polyhedron::find_shell_edge(const shell::Vert* v0, const shell::Vert* v1) const {
    for (auto& sedge : m_shell_edges)
        if ((sedge->verts[0] == v0 &&
             sedge->verts[1] == v1) ||
             (sedge->verts[1] == v0 &&
              sedge->verts[0] == v1))
            return sedge;

    return nullptr;
}


front::Face* Polyhedron::find_front_face(const Face* face) const {
    for (auto& fface : m_front_faces)
        if (fface->x == face)
            return fface;

    return nullptr;
}


std::vector<front::Edge*> Polyhedron::find_front_edge(const Vert* v0, const Vert* v1) const {
    std::vector<front::Edge*> res;
    for (auto& l_fedge : m_front_edges)
        if (((l_fedge->x->verts[0] == v0 &&
              l_fedge->x->verts[1] == v1) ||
              (l_fedge->x->verts[1] == v0 &&
               l_fedge->x->verts[0] == v1)))
            res.push_back(l_fedge);

    return res;
}


std::vector<front::Edge*> Polyhedron::find_front_edge(const Edge* edge) const {
    std::vector<front::Edge*> res;
    for (auto& l_fedge : m_front_edges)
        if (l_fedge->x == edge)
            res.push_back(l_fedge);

    return res;
}


front::Face* Polyhedron::add_to_front(const front::Face* fface, bool add_inner) {
    m_front_faces.push_back(const_cast<front::Face*>(fface));
    if (add_inner)
        m_inner_faces.push_back(fface->x);
    return const_cast<front::Face*>(fface);
}


front::Edge* Polyhedron::add_to_front(const front::Edge* fedge, bool add_inner) {
    m_front_edges.push_back(const_cast<front::Edge*>(fedge));
    if (add_inner)
        m_inner_edges.push_back(fedge->x);
    return const_cast<front::Edge*>(fedge);
}


front::Vert* Polyhedron::add_to_front(const front::Vert* fvert, bool add_inner) {
    m_front_verts.push_back(const_cast<front::Vert*>(fvert));
    if (add_inner)
        m_inner_verts.push_back(fvert->x);
    return const_cast<front::Vert*>(fvert);
}


front::Face* Polyhedron::add_to_front(const pmg::Face* face, bool add_inner) {
    auto new_fface = new front::Face(this, face);
    m_front_faces.push_back(new_fface);
    if (add_inner)
        m_inner_faces.push_back(new_fface->x);
    return new_fface;
}


front::Edge* Polyhedron::add_to_front(const pmg::Edge* edge, bool add_inner) {
    auto new_fedge = new front::Edge(this, edge);
    m_front_edges.push_back(new_fedge);
    if (add_inner)
        m_inner_edges.push_back(new_fedge->x);
    return new_fedge;
}

front::Vert* Polyhedron::add_to_front(const pmg::Vert* vert, bool add_inner) {
    auto new_fvert = new front::Vert(vert);
    m_front_verts.push_back(new_fvert);
    if (add_inner)
        m_inner_verts.push_back(new_fvert->x);
    return new_fvert;
}


void Polyhedron::remove_from_front(front::Face* fface) {
    //m_front_faces.erase(std::find(m_front_faces.begin(), m_front_faces.end(), fface));
    auto rem_iter = std::find(m_front_faces.begin(), m_front_faces.end(), fface);
    *rem_iter = m_front_faces.back();
    m_front_faces.pop_back();
}


void Polyhedron::remove_from_front(front::Edge* fedge) {
    //m_front_edges.erase(std::find(m_front_edges.begin(), m_front_edges.end(), fedge));
    auto rem_iter = std::find(m_front_edges.begin(), m_front_edges.end(), fedge);
    *rem_iter = m_front_edges.back();
    m_front_edges.pop_back();
}

void pmg::Polyhedron::remove_from_front(front::Vert* fvert) {
    auto rem_iter = std::find(m_front_verts.begin(), m_front_verts.end(), fvert);
    *rem_iter = m_front_verts.back();
    m_front_verts.pop_back();
}


bool Polyhedron::segmentIntersectMesh(const vec3& v0, const vec3& v1) const {
    for (auto& face : m_inner_faces) 
        if (spt::does_segment_intersect_triangle(
                v0, v1,
                face->edges[0]->verts[0]->pos(),
                face->edges[0]->verts[1]->pos(),
                face->find_vert_not(face->edges[0])->pos()))
            return true;

    return false;
}


bool Polyhedron::segmentIntersectFront(const vec3& v0, const vec3& v1) const {
    for (auto& fface : m_front_faces)
        if (spt::does_segment_intersect_triangle(
                v0, v1,
                fface->x->edges[0]->verts[0]->pos(),
                fface->x->edges[0]->verts[1]->pos(),
                fface->x->find_vert_not(fface->x->edges[0])->pos()))
            return true;

    return false;
}


bool Polyhedron::edge_intersect_front(const Vert* v0, const vec3& v1) const {
    for (auto& fface : m_front_faces) {
        if (fface->contains(v0))
            continue;

        if (spt::does_segment_intersect_triangle(
                v0->pos(), v1,
                fface->x->edges[0]->verts[0]->pos(),
                fface->x->edges[0]->verts[1]->pos(),
                fface->x->find_vert_not(fface->x->edges[0])->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::edge_intersect_front(const Vert* v0, const Vert* v1) const {
    for (auto& fface : m_front_faces) {
        if (fface->contains(v0)
            || fface->contains(v1))
            continue;

        if (spt::does_segment_intersect_triangle(
                v0->pos(), v1->pos(),
                fface->x->edges[0]->verts[0]->pos(),
                fface->x->edges[0]->verts[1]->pos(),
                fface->x->find_vert_not(fface->x->edges[0])->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::edge_intersect_any_face(const Edge* edge) const {
    for (auto& face : m_inner_faces) {
        if (face->contains(edge->verts[0])
            || face->contains(edge->verts[1]))
            continue;

        if (spt::does_segment_intersect_triangle(
                edge->verts[0]->pos(), edge->verts[1]->pos(),
                face->edges[0]->verts[0]->pos(),
                face->edges[0]->verts[1]->pos(),
                face->find_vert_not(face->edges[0])->pos()))
            return true;
    }

    return false;
}


bool XOR(bool b0, bool b1) {
    return (b0 || b1) && !(b0 && b1);
}


bool Polyhedron::will_edge_intersect_front(front::Edge* fedge) const {
    auto l_opp_verts = fedge->opp_verts();

    if (edge_intersect_front(std::get<0>(l_opp_verts), std::get<1>(l_opp_verts)))
        return true;

    for (auto& l_fedge : m_front_edges)
        if (l_fedge->x->contains(std::get<0>(l_opp_verts)) &&
            l_fedge->x->contains(std::get<1>(l_opp_verts)))
            return false;

    for (auto& l_fedge : m_front_edges) {
        Vert* vert_buf;
        std::array<bool, 2> contains{
            l_fedge->x->contains(std::get<0>(l_opp_verts)),
            l_fedge->x->contains(std::get<1>(l_opp_verts))
        };
        // TODO: change that doubtful checks
        if (contains[0]) {
            if (l_fedge->x->verts[0] == std::get<0>(l_opp_verts))
                vert_buf = l_fedge->x->verts[1];
            else
                vert_buf = l_fedge->x->verts[0];

            if (spt::distance_point_to_segment(vert_buf->pos(), std::get<0>(l_opp_verts)->pos(), std::get<1>(l_opp_verts)->pos()) < C_EDGES_INTERS_DIST * m_prefLen)
                return true;
        } else if (contains[1]) {
            if (l_fedge->x->verts[0] == std::get<1>(l_opp_verts))
                vert_buf = l_fedge->x->verts[1];
            else
                vert_buf = l_fedge->x->verts[0];

            if (spt::distance_point_to_segment(vert_buf->pos(), std::get<0>(l_opp_verts)->pos(), std::get<1>(l_opp_verts)->pos()) < C_EDGES_INTERS_DIST * m_prefLen)
                return true;
        } else {
            if (spt::segments_distance(
                    std::get<0>(l_opp_verts)->pos(), std::get<1>(l_opp_verts)->pos(),
                    l_fedge->x->verts[0]->pos(), l_fedge->x->verts[1]->pos()) < C_EDGES_INTERS_DIST * m_prefLen)
                return true;
        }
    }

    return false;
}


bool Polyhedron::any_edge_intersect_face(const Vert* v0, const Vert* v1, const vec3& v2) const {
    for (auto& l_fedge : m_front_edges) {
        std::array<bool, 2> contains{
            l_fedge->x->contains(v0),
            l_fedge->x->contains(v1)
        };

        if (contains[0] || contains[1])
            continue;

        auto cur_verts = l_fedge->x->verts;
        if (spt::does_segment_intersect_triangle(
                cur_verts[0]->pos(), cur_verts[1]->pos(),
                v0->pos(), v1->pos(), v2))
            return true;
    }

    return false;
}


bool Polyhedron::any_edge_intersect_face(const Vert* v0, const Vert* v1, const Vert* v2) const {
    for (auto& l_fedge : m_front_edges) {
        if (l_fedge->x->contains(v0)
            || l_fedge->x->contains(v1)
            || l_fedge->x->contains(v2))
            continue;

        if (spt::does_segment_intersect_triangle(
                l_fedge->x->verts[0]->pos(), l_fedge->x->verts[1]->pos(),
                v0->pos(), v1->pos(), v2->pos()))
            return true;
    }

    return false;
}


bool Polyhedron::will_any_edge_intersect_faces(front::Edge* fedge) const {
    auto l_opp_verts = fedge->opp_verts();

    if (any_edge_intersect_face(fedge->x->verts[0], std::get<0>(l_opp_verts), std::get<1>(l_opp_verts)) ||
        any_edge_intersect_face(fedge->x->verts[1], std::get<0>(l_opp_verts), std::get<1>(l_opp_verts)))
        return true;

    return false;
}


bool Polyhedron::will_any_vert_inside_tetr(front::Edge* fedge) const {
    auto l_opp_verts = fedge->opp_verts();

    std::array<vec3, 4> points;
    points[0] = std::get<0>(l_opp_verts)->pos();
    points[1] = std::get<1>(l_opp_verts)->pos();
    points[2] = fedge->x->verts[0]->pos();
    points[3] = fedge->x->verts[1]->pos();
    for (auto& vert : m_inner_verts) {
        if (vert != std::get<0>(l_opp_verts) &&
            vert != std::get<1>(l_opp_verts) &&
            vert != fedge->x->verts[0] &&
            vert != fedge->x->verts[1] &&
            spt::is_point_in_tetrahedron(vert->pos(), points[0], points[1], points[2], points[3]))
            return true;
    }

    return false;
}


bool Polyhedron::will_front_split(front::Edge* fedge, front::Edge* opp_fedge) const {
    auto l_opp_fedge = opp_fedge;
    if (!l_opp_fedge)
        l_opp_fedge = fedge->find_opp_edge();
    if (!l_opp_fedge)
        return false;

    auto opp_fedge_opp_verts = l_opp_fedge->opp_verts();

    if (std::get<0>(opp_fedge_opp_verts) == fedge->x->verts[0] ||
        std::get<0>(opp_fedge_opp_verts) == fedge->x->verts[1] ||
        std::get<1>(opp_fedge_opp_verts) == fedge->x->verts[0] ||
        std::get<1>(opp_fedge_opp_verts) == fedge->x->verts[1])
        return false;

    return true;
}


bool Polyhedron::will_front_collapse(front::Edge* fedge, front::Edge* opp_fedge) const {
    auto l_opp_fedge = opp_fedge;
    if (!l_opp_fedge)
        l_opp_fedge = fedge->find_opp_edge();
    if (!l_opp_fedge)
        return false;

    auto opp_ffaces = l_opp_fedge->adj_ffaces();

    if (!fedge->x->contains(opp_ffaces.first->x->find_vert_not(l_opp_fedge->x)) ||
        !fedge->x->contains(opp_ffaces.second->x->find_vert_not(l_opp_fedge->x)))
        return false;

    return true;
}


bool Polyhedron::will_parallel_faces(front::Edge* fedge) const {
    auto adj_ffaces = fedge->adj_ffaces();

    std::array<Vert*, 2> l_opp_verts;
    l_opp_verts[0] = std::get<0>(adj_ffaces)->x->find_vert_not(fedge->x);
    l_opp_verts[1] = std::get<1>(adj_ffaces)->x->find_vert_not(fedge->x);

    std::array<vec3, 2> plane0;
    plane0[0] = l_opp_verts[0]->pos() - fedge->x->verts[0]->pos();
    plane0[1] = l_opp_verts[1]->pos() - fedge->x->verts[0]->pos();
    vec3 normal0 = spt::cross(plane0[0], plane0[1]).normalize();
    std::array<vec3, 2> plane1;
    plane1[0] = l_opp_verts[0]->pos() - fedge->x->verts[1]->pos();
    plane1[1] = l_opp_verts[1]->pos() - fedge->x->verts[1]->pos();
    vec3 normal1 = spt::cross(plane1[0], plane1[1]).normalize();

    for (auto& fface : m_front_faces) {
        std::array<Edge*, 2> inters_reses;
        if ((fface != std::get<0>(adj_ffaces)) &&
            (fface != std::get<1>(adj_ffaces))) {
            inters_reses[0] = relations::adjacent_by_edge(fface->x, std::get<0>(adj_ffaces)->x);
            inters_reses[1] = relations::adjacent_by_edge(fface->x, std::get<1>(adj_ffaces)->x);
            if (!XOR(static_cast<bool>(inters_reses[0]), static_cast<bool>(inters_reses[1])))
                continue;

            std::array<Vert*, 2> fface_to_verts;
            fface_to_verts[0] = fface->x->edges[0]->verts[0];
            fface_to_verts[1] = fface->x->edges[0]->verts[1];
            Vert* fface_from_vert = fface->x->find_vert_not(fface->x->edges[0]);

            std::array<vec3, 2> f_plane;
            f_plane[0] = fface_to_verts[0]->pos() - fface_from_vert->pos();
            f_plane[1] = fface_to_verts[1]->pos() - fface_from_vert->pos();
            vec3 f_normal = spt::cross(f_plane[0], f_plane[1]).normalize();

            if (std::abs(std::abs(spt::dot(f_normal, normal0)) - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6) ||
                std::abs(std::abs(spt::dot(f_normal, normal1)) - static_cast<real_t>(1.0)) < static_cast<real_t>(1e-6)) {
                std::size_t i = inters_reses[0] ? 0 : 1;

                std::array<vec3, 2> border_verts;
                border_verts[0] = inters_reses[i]->verts[0]->pos();
                border_verts[1] = inters_reses[i]->verts[1]->pos();

                vec3 main_face_3rd_vert;
                if (inters_reses[i]->contains(l_opp_verts[0]))
                    main_face_3rd_vert = l_opp_verts[1]->pos();
                else
                    main_face_3rd_vert = l_opp_verts[0]->pos();

                vec3 cur_face_3rd_vert = fface->x->find_vert_not(inters_reses[i])->pos();

                vec3 main_face_cross = spt::cross(main_face_3rd_vert - border_verts[0], main_face_3rd_vert - border_verts[1]);
                vec3 cur_face_cross = spt::cross(cur_face_3rd_vert - border_verts[0], cur_face_3rd_vert - border_verts[1]);
                if (spt::dot(main_face_cross, cur_face_cross) > static_cast<real_t>(0.0))
                    return true;
            }
        }
    }

    return false;
}


bool Polyhedron::does_front_intersect_sphere(const vec3& center, real_t radius) const {
    for (auto& fface : m_front_faces) {
        std::array<vec3, 3> triangle;
        triangle[0] = fface->x->edges[0]->verts[0]->pos();
        triangle[1] = fface->x->edges[0]->verts[1]->pos();
        triangle[2] = fface->x->find_vert_not(fface->x->edges[0])->pos();
        if (spt::does_triangle_intersect_sphere(triangle[0], triangle[1], triangle[2], center, radius))
            return true;
    }

    return false;
}


std::pair<real_t, real_t> Polyhedron::min_max_edges_lengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3) {
    auto min_max = min_max_edges_sqr_lengths(p0, p1, p2, p3);
    min_max.first = std::sqrt(min_max.first);
    min_max.second = std::sqrt(min_max.second);
    return min_max;
}


std::pair<real_t, real_t> Polyhedron::min_max_edges_sqr_lengths(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3) {
    std::array<real_t, 6> sqr_magns;
    sqr_magns[0] = (p1 - p0).sqr_magnitude();
    sqr_magns[1] = (p2 - p0).sqr_magnitude();
    sqr_magns[2] = (p3 - p0).sqr_magnitude();
    sqr_magns[3] = (p2 - p1).sqr_magnitude();
    sqr_magns[4] = (p3 - p1).sqr_magnitude();
    sqr_magns[5] = (p3 - p2).sqr_magnitude();
    return std::minmax({ sqr_magns[0], sqr_magns[1], sqr_magns[2], sqr_magns[3], sqr_magns[4], sqr_magns[5] });
}


real_t Polyhedron::tetr_simple_quality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3) {
    auto sqr_min_max = min_max_edges_sqr_lengths(p0, p1, p2, p3);
    return std::sqrt(sqr_min_max.first / sqr_min_max.second);
}


real_t Polyhedron::tetr_simple_sqr_quality(const vec3& p0, const vec3& p1, const vec3& p2, const vec3& p3) {
    auto sqr_min_max = min_max_edges_sqr_lengths(p0, p1, p2, p3);
    return sqr_min_max.first / sqr_min_max.second;
}


front::Edge* Polyhedron::current_front_edge(real_t max_compl) const {
    real_t cur_max_compl = static_cast<real_t>(0);
    front::Edge* max_fedge = nullptr;
    for (auto& l_fedge : m_front_edges) {
        real_t cur_compl = l_fedge->complexity();
        if (cur_compl > cur_max_compl &&
            cur_compl < max_compl) {
            cur_max_compl = cur_compl;
            max_fedge = l_fedge;
        }
    }

    return max_fedge;
}


bool Polyhedron::exhaust_without_new_vert_priority_predicate(front::Edge* curFEdge) {
    if (curFEdge->angle() < degToRad(70))
        return true;

    auto l_opp_verts = curFEdge->opp_verts();
    if (curFEdge->find_opp_edge() ||
        (curFEdge->angle() > degToRad(70) &&
         curFEdge->angle() < degToRad(100) &&
         (std::get<1>(l_opp_verts)->pos() - std::get<0>(l_opp_verts)->pos()).sqr_magnitude() <= m_prefLen * m_prefLen))
        return true;

    if (will_front_split(curFEdge))
        return true;

    return false;
}


bool Polyhedron::exhaust_with_new_vert_priority_predicate(front::Edge* current_fedge) {
    if (current_fedge->angle() > degToRad(120))
        return true;

    return false;
}


Polyhedron::exhaust_type Polyhedron::exhaustion_type_quality_priority(
    front::Edge* current_fedge,
    front::Face*& out_with_nv_fface, vec3*& out_with_nv_new_vert_pos) {
    if (will_front_split(current_fedge))
        return exhaust_type::without_new_vert;

    if (will_parallel_faces(current_fedge) || // TODO: replace later with check if front is close to potential element
        will_any_edge_intersect_faces(current_fedge) ||
        will_any_vert_inside_tetr(current_fedge))
        return exhaust_type::with_new_vert;

    auto l_opp_verts = current_fedge->opp_verts();
    real_t without_nv_quality = tetr_simple_sqr_quality(
        current_fedge->x->verts[0]->pos(),
        current_fedge->x->verts[1]->pos(),
        std::get<0>(l_opp_verts)->pos(),
        std::get<1>(l_opp_verts)->pos());

    front::Face* fface = choose_fface_for_exhaustion_with_new_vert(current_fedge);
    vec3 new_vert_pos;
    if (!try_compute_new_vert_pos(fface, new_vert_pos))
        return exhaust_type::dont_exhaust;

    real_t with_nv_quality = tetr_simple_sqr_quality(
        fface->x->edges[0]->verts[0]->pos(),
        fface->x->edges[0]->verts[1]->pos(),
        fface->x->find_vert_not(fface->x->edges[0])->pos(),
        new_vert_pos);

    if (without_nv_quality > with_nv_quality)
        return exhaust_type::without_new_vert;

    out_with_nv_fface = fface;
    out_with_nv_new_vert_pos = new vec3(new_vert_pos);
    return exhaust_type::with_new_vert;
}


vec3 Polyhedron::compute_normal_in_tetr(const front::Face* fface, const vec3& opp_vert_pos) const {
    return compute_normal_in_tetr(fface->x->edges[0]->verts[0]->pos(),
                                  fface->x->edges[0]->verts[1]->pos(),
                                  fface->x->find_vert_not(fface->x->edges[0])->pos(),
                                  opp_vert_pos);
}


vec3 Polyhedron::compute_normal_in_tetr(const front::Face* fface, const pmg::Edge* one_of_remaining_edges) const {
    vec3 opp_pos = fface->x->contains(one_of_remaining_edges->verts[0]) ?
        one_of_remaining_edges->verts[1]->pos() :
        one_of_remaining_edges->verts[0]->pos();

    return compute_normal_in_tetr(fface, opp_pos);
}


vec3 Polyhedron::compute_normal_in_tetr(const vec3& fface_pos0, const vec3& fface_pos1, const vec3& fface_pos2, const vec3& opp_vert_pos) const {
    vec3 normal = spt::cross(
        fface_pos0 - fface_pos2,
        fface_pos1 - fface_pos2).normalize();
    if (spt::dot(normal, opp_vert_pos - fface_pos2) > static_cast<real_t>(0.0))
        normal *= static_cast<real_t>(-1.0);

    return normal;
}


void Polyhedron::set_front_edges_in_front_split(const front::Edge* fedge, std::array<front::Edge*, 2> new_opp_fedges, std::array<front::Face*, 2> new_ffaces, pair_ff opp_ffaces) const {
    pmg::Edge* opp_edge = new_opp_fedges[0]->x;
    std::array<vec3, 2> opp_verts_poses;
    opp_verts_poses[0] = opp_edge->verts[0]->pos();
    opp_verts_poses[1] = opp_edge->verts[1]->pos();
    vec3 main_vert0_pos = new_ffaces[0]->x->contains(fedge->x->verts[0]) ?
        fedge->x->verts[0]->pos() : fedge->x->verts[1]->pos();
    vec3 main_vert1_pos = new_ffaces[1]->x->contains(fedge->x->verts[0]) ?
        fedge->x->verts[0]->pos() : fedge->x->verts[1]->pos();
    vec3 main_vert0_proj = spt::project(main_vert0_pos, opp_verts_poses[0], opp_verts_poses[1]);
    vec3 main_vert1_proj = spt::project(main_vert1_pos, opp_verts_poses[0], opp_verts_poses[1]);
    vec3 main_vec0 = main_vert0_pos - main_vert0_proj;
    vec3 main_vec1 = main_vert1_pos - main_vert1_proj;

    vec3 adj_opp_pos0 = opp_ffaces.first->x->find_vert_not(opp_edge)->pos();
    vec3 adj_opp_pos1 = opp_ffaces.second->x->find_vert_not(opp_edge)->pos();
    vec3 adj_opp_proj0 = spt::project(adj_opp_pos0, opp_verts_poses[0], opp_verts_poses[1]);
    vec3 adj_opp_proj1 = spt::project(adj_opp_pos1, opp_verts_poses[0], opp_verts_poses[1]);
    vec3 adj_vec0 = adj_opp_pos0 - adj_opp_proj0;
    vec3 adj_vec1 = adj_opp_pos1 - adj_opp_proj1;

    std::array<std::array<real_t, 2>, 2> coses;
    coses[0][0] = spt::cos(main_vec0, adj_vec0);
    coses[0][1] = spt::cos(main_vec0, adj_vec1);
    coses[1][0] = spt::cos(main_vec1, adj_vec0);
    coses[1][1] = spt::cos(main_vec1, adj_vec1);

    if (coses[0][0] > coses[0][1] && coses[1][1] > coses[1][0]) {
        new_opp_fedges[0]->fill_adj_ffaces(new_ffaces[0], opp_ffaces.first);
        opp_ffaces.first->add_front_edge(new_opp_fedges[0]);

        new_opp_fedges[1]->fill_adj_ffaces(new_ffaces[1], opp_ffaces.second);
        opp_ffaces.second->add_front_edge(new_opp_fedges[1]);
        return;
    }

    if (coses[0][1] >= coses[0][0] && coses[1][0] >= coses[1][1]) {
        new_opp_fedges[0]->fill_adj_ffaces(new_ffaces[0], opp_ffaces.second);
        opp_ffaces.second->add_front_edge(new_opp_fedges[0]);

        new_opp_fedges[1]->fill_adj_ffaces(new_ffaces[1], opp_ffaces.first);
        opp_ffaces.first->add_front_edge(new_opp_fedges[1]);
        return;
    }

    std::size_t best_ofi;
    if (coses[0][0] > coses[0][1] && coses[1][0] > coses[1][1])
        best_ofi = 0;
    else
        // When coses[0][1] >= coses[0][0] && coses[1][1] >= coses[1][0]
        best_ofi = 1;

    std::size_t worst_ofi = best_ofi == 0 ? 1 : 0;

    std::array<front::Face*, 2> l_opp_ffaces = { opp_ffaces.first, opp_ffaces.second };
    if (coses[0][best_ofi] > coses[1][best_ofi]) {
        new_opp_fedges[0]->fill_adj_ffaces(new_ffaces[0], l_opp_ffaces[best_ofi]);
        l_opp_ffaces[best_ofi]->add_front_edge(new_opp_fedges[0]);

        new_opp_fedges[1]->fill_adj_ffaces(new_ffaces[1], l_opp_ffaces[worst_ofi]);
        l_opp_ffaces[worst_ofi]->add_front_edge(new_opp_fedges[1]);
    } else {
        new_opp_fedges[0]->fill_adj_ffaces(new_ffaces[0], l_opp_ffaces[worst_ofi]);
        l_opp_ffaces[worst_ofi]->add_front_edge(new_opp_fedges[0]);

        new_opp_fedges[1]->fill_adj_ffaces(new_ffaces[1], l_opp_ffaces[best_ofi]);
        l_opp_ffaces[best_ofi]->add_front_edge(new_opp_fedges[1]);
    }

    return;
}


void Polyhedron::exhaust_front_collapse(front::Edge* fedge, front::Edge* opp_fedge) {
    auto fedge_adj_faces = fedge->adj_ffaces();
    auto opp_ffaces = opp_fedge->adj_ffaces();

    m_inner_tetrs.push_back(new Tetr(
        fedge->x->verts[0],
        fedge->x->verts[1],
        opp_fedge->x->verts[0],
        opp_fedge->x->verts[1]));

    auto new_tetr = new Tetr(fedge->x->verts[0],
                             fedge->x->verts[1],
                             opp_fedge->x->verts[0],
                             opp_fedge->x->verts[1]);
    m_inner_tetrs.push_back(new_tetr);
    
    std::vector<front::Edge*> fedges_to_erase;
    fedges_to_erase.reserve(4);
    for (auto& l_fedge : opp_ffaces.first->front_edges)
        if (l_fedge->x->contains(fedge->x->verts[0]) || l_fedge->x->contains(fedge->x->verts[1]))
            fedges_to_erase.push_back(l_fedge);

    for (auto& l_fedge : opp_ffaces.second->front_edges)
        if (l_fedge->x->contains(fedge->x->verts[0]) || l_fedge->x->contains(fedge->x->verts[1]))
            fedges_to_erase.push_back(l_fedge);


    remove_from_front(fedge_adj_faces.first);
    delete fedge_adj_faces.first;
    remove_from_front(fedge_adj_faces.second);
    delete fedge_adj_faces.second;
    remove_from_front(opp_ffaces.first);
    delete opp_ffaces.first;
    remove_from_front(opp_ffaces.second);
    delete opp_ffaces.second;
    remove_from_front(fedge);
    delete fedge;
    remove_from_front(opp_fedge);
    delete opp_fedge;
    for (auto& l_fedge : fedges_to_erase) {
        remove_from_front(l_fedge);
        delete l_fedge;
    }

#ifdef DEBUG
    debug();
    std::cout << std::endl << "exhaust_front_collapse";
#endif // DEBUG
}


void Polyhedron::exhaust_front_split(front::Edge* fedge, front::Edge* opp_fedge) {
    // HACK: this volatile helps to avoid computational error
    volatile auto adj_ffaces = fedge->adj_ffaces();

    std::array<front::Vert*, 2> l_opp_fverts;
    l_opp_fverts[0] = adj_ffaces.first->find_front_vert_not(fedge);
    l_opp_fverts[1] = adj_ffaces.second->find_front_vert_not(fedge);
    
    auto opp_ffaces = opp_fedge->adj_ffaces();
    opp_ffaces.first->remove_front_edge(opp_fedge);
    opp_ffaces.second->remove_front_edge(opp_fedge);
    remove_from_front(opp_fedge);

    std::array<std::array<front::Vert*, 2>, 2> new_opp_fverts;
    new_opp_fverts[0][0] = add_to_front(new front::Vert(opp_fedge->x->verts[0]), false);
    new_opp_fverts[0][1] = add_to_front(new front::Vert(opp_fedge->x->verts[1]), false);
    new_opp_fverts[1][0] = add_to_front(new front::Vert(opp_fedge->x->verts[0]), false);
    new_opp_fverts[1][1] = add_to_front(new front::Vert(opp_fedge->x->verts[1]), false);

    std::array<front::Edge*, 2> new_opp_fedges;
    new_opp_fedges[0] = add_to_front(new front::Edge(this, new_opp_fverts[0][0], new_opp_fverts[0][1]), false);
    new_opp_fedges[1] = add_to_front(new front::Edge(this, new_opp_fverts[1][0], new_opp_fverts[1][1]), false);

    ////////

    std::array<front::Edge*, 3> new_fface_fedges;
    new_fface_fedges[0] = nullptr;
    new_fface_fedges[1] = nullptr;

    std::array<front::Face*, 2> new_ffaces;

    new_fface_fedges[2] = new_opp_fedges[0];
    for (auto& l_fedge : adj_ffaces.first->front_edges) {
        if (l_fedge->contains(fedge->front_verts[0]) &&
            (l_fedge->contains(opp_fedge->front_verts[0]) ||
             l_fedge->contains(opp_fedge->front_verts[1]))) {
            new_fface_fedges[0] = l_fedge;
            l_fedge->refresh_angle_data();
            l_fedge->remove_adj_fface(adj_ffaces.first);
            l_fedge->remove_adj_fface(adj_ffaces.second);
            break;
        }
    }
    for (auto& l_fedge : adj_ffaces.second->front_edges) {
        if (l_fedge->contains(fedge->front_verts[0]) &&
            (l_fedge->contains(opp_fedge->front_verts[0]) ||
             l_fedge->contains(opp_fedge->front_verts[1]))) {
            new_fface_fedges[1] = l_fedge;
            l_fedge->refresh_angle_data();
            l_fedge->remove_adj_fface(adj_ffaces.first);
            l_fedge->remove_adj_fface(adj_ffaces.second);
            break;
        }
    }
    new_ffaces[0] = add_to_front(new front::Face(this, new_fface_fedges[0],
                                                       new_fface_fedges[1],
                                                       new_fface_fedges[2]));
    new_ffaces[0]->normal = compute_normal_in_tetr(new_ffaces[0], fedge->x);


    new_fface_fedges[2] = new_opp_fedges[1];
    for (auto& l_fedge : adj_ffaces.first->front_edges) {
        if (l_fedge->contains(fedge->front_verts[1]) &&
            (l_fedge->contains(opp_fedge->front_verts[0]) ||
             l_fedge->contains(opp_fedge->front_verts[1]))) {
            new_fface_fedges[0] = l_fedge;
            l_fedge->refresh_angle_data();
            l_fedge->remove_adj_fface(adj_ffaces.first);
            l_fedge->remove_adj_fface(adj_ffaces.second);
            break;
        }
    }
    for (auto& l_fedge : adj_ffaces.second->front_edges) {
        if (l_fedge->contains(fedge->front_verts[1]) &&
            (l_fedge->contains(opp_fedge->front_verts[0]) ||
             l_fedge->contains(opp_fedge->front_verts[1]))) {
            new_fface_fedges[1] = l_fedge;
            l_fedge->refresh_angle_data();
            l_fedge->remove_adj_fface(adj_ffaces.first);
            l_fedge->remove_adj_fface(adj_ffaces.second);
            break;
        }
    }
    new_ffaces[1] = add_to_front(new front::Face(this, new_fface_fedges[0],
                                                       new_fface_fedges[1],
                                                       new_fface_fedges[2]));
    new_ffaces[1]->normal = compute_normal_in_tetr(new_ffaces[1], fedge->x);


    set_front_edges_in_front_split(fedge, { new_opp_fedges[0], new_opp_fedges[1] }, { new_ffaces[0], new_ffaces[1] }, opp_ffaces);
    
    m_inner_tetrs.push_back(new Tetr(fedge->front_verts[0]->x,
                                     fedge->front_verts[1]->x,
                                     l_opp_fverts[0]->x,
                                     l_opp_fverts[1]->x));
    

    delete opp_fedge;
    remove_from_front(adj_ffaces.first);
    delete adj_ffaces.first;
    remove_from_front(adj_ffaces.second);
    delete adj_ffaces.second;
    remove_from_front(fedge);
    delete fedge;

#ifdef DEBUG
    debug();
    std::cout << std::endl << "exhaust_front_split";
#endif // DEBUG
}


void Polyhedron::exhaust_without_new_vert_opp_edge_exists(front::Edge* fedge, front::Edge* opp_fedge) {
#ifdef DEBUG
    debug();
    opp_fedge = fedge->find_opp_edge();
    opp_fedge->clear_adj_ffaces();
#endif // DEBUG
    auto opp_ffaces = opp_fedge->adj_ffaces();

#ifdef DEBUG
    debug();
    fedge->clear_adj_ffaces();
#endif // DEBUG
    
    auto fedge_adj_faces = fedge->adj_ffaces();
    std::array<front::Face*, 3> main_ffaces;
    main_ffaces[0] = std::get<0>(fedge_adj_faces);
    main_ffaces[1] = std::get<1>(fedge_adj_faces);

    front::Vert* main_fvert;
    if (std::get<0>(opp_ffaces)->contains(fedge->front_verts[0])) {
        main_ffaces[2] = std::get<0>(opp_ffaces);
        main_fvert = fedge->front_verts[0];
    } else if (std::get<0>(opp_ffaces)->x->contains(fedge->x->verts[1])) { //
        main_ffaces[2] = std::get<0>(opp_ffaces);
        main_fvert = fedge->front_verts[1];
    } else if (std::get<1>(opp_ffaces)->x->contains(fedge->x->verts[0])) {
        main_ffaces[2] = std::get<1>(opp_ffaces);
        main_fvert = fedge->front_verts[0];
    } else {
        main_ffaces[2] = std::get<1>(opp_ffaces);
        main_fvert = fedge->front_verts[1];
    }

#ifdef DEBUG
    debug();
    std::array<front::Edge*, 3> fedges_to_erase_deb;
    std::size_t idx_deb = 0;
    for (auto& l_fface : main_ffaces)
        for (auto& l_fedge : l_fface->front_edges)
            if (std::find(fedges_to_erase_deb.begin(), fedges_to_erase_deb.begin() + idx_deb, l_fedge) 
                == fedges_to_erase_deb.begin() + idx_deb &&
                l_fedge->contains(main_fvert))
                fedges_to_erase_deb[idx_deb++] = l_fedge;
#endif // DEBUG

    std::array<front::Edge*, 3> new_tetr_fedges;
    for (std::size_t i = 0; i < 3; i++) {
        new_tetr_fedges[i] = main_ffaces[i]->find_front_edge_not(main_fvert);
        new_tetr_fedges[i]->refresh_angle_data();
        new_tetr_fedges[i]->remove_adj_fface(main_ffaces[i]);
    }

    auto new_fface = add_to_front(new front::Face(this, new_tetr_fedges[0],
                                                        new_tetr_fedges[1],
                                                        new_tetr_fedges[2]));
    new_fface->normal = compute_normal_in_tetr(new_fface, main_fvert->x->pos());
    
    m_inner_tetrs.push_back(new Tetr(fedge->front_verts[0]->x,
                                     fedge->front_verts[1]->x,
                                     opp_fedge->front_verts[0]->x,
                                     opp_fedge->front_verts[1]->x));
        
    std::array<front::Edge*, 3> fedges_to_erase;
    std::size_t idx = 0;
    for (auto& l_fface : main_ffaces)
        for (auto& l_fedge : l_fface->front_edges)
            if (std::find(fedges_to_erase.begin(), fedges_to_erase.begin() + idx, l_fedge) == fedges_to_erase.begin() + idx &&
                l_fedge->contains(main_fvert))
                fedges_to_erase[idx++] = l_fedge;

    for (auto& l_fedge : fedges_to_erase) {
        remove_from_front(l_fedge);
        delete l_fedge;
    }

    for (auto& l_fface : main_ffaces) {
        remove_from_front(l_fface);
        delete l_fface;
    }

#ifdef DEBUG
    debug();
    std::cout << std::endl << "exhaust_without_new_vert_opp_edge_exists";
#endif // DEBUG
}


void Polyhedron::exhaust_without_new_vert_opp_edge_doesnt_exist(front::Edge* fedge) {
    // HACK: this volatile helps to avoid computational error.
    volatile auto adj_ffaces = fedge->adj_ffaces();

    std::array<front::Vert*, 2> l_opp_fverts;
    l_opp_fverts[0] = adj_ffaces.first->find_front_vert_not(fedge);
    l_opp_fverts[1] = adj_ffaces.second->find_front_vert_not(fedge);

    front::Edge* opp_fedge = add_to_front(new front::Edge(this, l_opp_fverts[0], l_opp_fverts[1]));

    std::array<front::Edge*, 3> new_fface_fedges;
    new_fface_fedges[0] = nullptr;
    new_fface_fedges[1] = nullptr;
    new_fface_fedges[2] = opp_fedge;

    for (auto& l_fedge : adj_ffaces.first->front_edges) {
        if (l_fedge->contains(fedge->front_verts[0]) &&
            (l_fedge->contains(opp_fedge->front_verts[0]) ||
             l_fedge->contains(opp_fedge->front_verts[1]))) {
            new_fface_fedges[0] = l_fedge;
            l_fedge->refresh_angle_data();
            l_fedge->remove_adj_fface(adj_ffaces.first);
//            l_fedge->remove_adj_fface(adj_ffaces.second);
            break;
        }
    }
    for (auto& l_fedge : adj_ffaces.second->front_edges) {
        if (l_fedge->contains(fedge->front_verts[0]) &&
            (l_fedge->contains(opp_fedge->front_verts[0]) ||
             l_fedge->contains(opp_fedge->front_verts[1]))) {
            new_fface_fedges[1] = l_fedge;
            l_fedge->refresh_angle_data();
//            l_fedge->remove_adj_fface(adj_ffaces.first);
            l_fedge->remove_adj_fface(adj_ffaces.second);
            break;
        }
    }
    auto new_fface = add_to_front(new front::Face(this, new_fface_fedges[0],
                                                        new_fface_fedges[1],
                                                        new_fface_fedges[2]));
    new_fface->normal = compute_normal_in_tetr(new_fface, fedge->x);

    for (auto& l_fedge : adj_ffaces.first->front_edges) {
        if (l_fedge->contains(fedge->front_verts[1]) &&
            (l_fedge->contains(opp_fedge->front_verts[0]) ||
             l_fedge->contains(opp_fedge->front_verts[1]))) {
            new_fface_fedges[0] = l_fedge;
            l_fedge->refresh_angle_data();
            l_fedge->remove_adj_fface(adj_ffaces.first);
            l_fedge->remove_adj_fface(adj_ffaces.second);
            break;
        }
    }
    for (auto& l_fedge : adj_ffaces.second->front_edges) {
        if (l_fedge->contains(fedge->front_verts[1]) &&
            (l_fedge->contains(opp_fedge->front_verts[0]) ||
             l_fedge->contains(opp_fedge->front_verts[1]))) {
            new_fface_fedges[1] = l_fedge;
            l_fedge->refresh_angle_data();
            l_fedge->remove_adj_fface(adj_ffaces.first);
            l_fedge->remove_adj_fface(adj_ffaces.second);
            break;
        }
    }
    new_fface = add_to_front(new front::Face(this, new_fface_fedges[0],
                                                   new_fface_fedges[1],
                                                   new_fface_fedges[2]));
    new_fface->normal = compute_normal_in_tetr(new_fface, fedge->x);
    
    m_inner_tetrs.push_back(new Tetr(fedge->front_verts[0]->x,
                                     fedge->front_verts[1]->x,
                                     l_opp_fverts[0]->x,
                                     l_opp_fverts[1]->x));
    
    remove_from_front(adj_ffaces.first);
    delete adj_ffaces.first;
    remove_from_front(adj_ffaces.second);
    delete adj_ffaces.second;
    remove_from_front(fedge);
    delete fedge;

#ifdef DEBUG
    debug();
    std::cout << std::endl << "exhaust_without_new_vert_opp_edge_doesnt_exist";
#endif // DEBUG
}


void Polyhedron::exhaust_without_new_vert(front::Edge* fedge, bool opp_edge_existence, front::Edge* opp_fedge) {
    front::Edge* l_opp_fedge = nullptr;
    if (opp_edge_existence && opp_fedge)
        l_opp_fedge = opp_fedge;
    else if (opp_edge_existence)
        l_opp_fedge = fedge->find_opp_edge();

    if (l_opp_fedge) {
        if (will_front_collapse(fedge, l_opp_fedge))
            exhaust_front_collapse(fedge, l_opp_fedge);

        else if (will_front_split(fedge, l_opp_fedge))
            exhaust_front_split(fedge, l_opp_fedge);

        else
            exhaust_without_new_vert_opp_edge_exists(fedge, l_opp_fedge);
    } else {
        exhaust_without_new_vert_opp_edge_doesnt_exist(fedge);
    }
}


bool Polyhedron::try_compute_new_vert_pos_type3(front::Face* fface, vec3& out_pos) {
    std::array<front::Edge*, 2> main_fedges;
    main_fedges[0] = fface->front_edges[0];
    main_fedges[1] = fface->front_edges[1];
    std::array<Edge*, 2> main_edges;
    main_edges[0] = main_fedges[0]->x;
    main_edges[1] = main_fedges[1]->x;

    auto main_vert = main_edges[0]->verts[0] == main_edges[1]->verts[0] ?
        main_edges[0]->verts[0] :
        (main_edges[0]->verts[0] == main_edges[1]->verts[1] ?
         main_edges[0]->verts[0] :
         main_edges[0]->verts[1]);
    auto third_edge = fface->x->find_edge_not(main_vert);
    auto third_f_edge = fface->find_front_edge(third_edge);

    auto v0 = main_edges[0]->verts[0] == main_vert ? main_edges[0]->verts[1] : main_edges[0]->verts[0];
    auto v1 = main_edges[1]->verts[0] == main_vert ? main_edges[1]->verts[1] : main_edges[1]->verts[0];
    auto v2 = main_vert;

    auto adj_ffaces = main_fedges[0]->adj_ffaces();
    auto fn0 = std::get<0>(adj_ffaces) == fface ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    adj_ffaces = main_fedges[1]->adj_ffaces();
    auto fn1 = std::get<0>(adj_ffaces) == fface ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    adj_ffaces = third_f_edge->adj_ffaces();
    auto fn2 = std::get<0>(adj_ffaces) == fface ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);

    vec3 v0_pos = v0->pos();
    vec3 v1_pos = v1->pos();
    vec3 v2_pos = v2->pos();

    vec3 n_m = fface->normal;
    vec3 n_n0 = fn0->normal;
    vec3 n_n1 = fn1->normal;
    vec3 n_n2 = fn2->normal;

    vec3 e_mn0 = (n_m + n_n0).normalize();
    vec3 e_mn1 = (n_m + n_n1).normalize();
    vec3 e_mn2 = n_m + n_n2;

    vec3 np_mn0 = spt::cross(v2_pos - v0_pos, e_mn0);
    vec3 np_mn1 = spt::cross(v2_pos - v1_pos, e_mn1);

    vec3 e = spt::cross(np_mn0, np_mn1);

    vec3 new_pos = spt::line_intersect_plane(v2_pos, e, v0_pos, v1_pos, v0_pos + e_mn2);

    std::array<real_t, 4> sum_magns;
    std::size_t i = 0;
    for (auto& fface : { fn0, fn1, fn2, fface })
        sum_magns[i++] = fface->x->edges[0]->magnitude()
        + fface->x->edges[1]->magnitude()
        + fface->x->edges[2]->magnitude();
    real_t av_magn = (sum_magns[0] + sum_magns[1] + sum_magns[2] + sum_magns[3]) / static_cast<real_t>(12.0);
    if (edge_intersect_front(v0, new_pos) ||
        edge_intersect_front(v1, new_pos) ||
        edge_intersect_front(v2, new_pos) ||
        does_front_intersect_sphere(new_pos, C_MIN_DIS * av_magn) ||
        any_edge_intersect_face(v0, v1, new_pos) ||
        any_edge_intersect_face(v0, v2, new_pos) ||
        any_edge_intersect_face(v1, v2, new_pos))
        return false;

//    std::cout << std::endl << "Type3";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::try_compute_new_vert_pos_type2(front::Face* fface, vec3& out_pos, std::size_t small_angle_idx0, std::size_t small_angle_idx1) {
    std::array<front::Edge*, 2> main_fedges;
    main_fedges[0] = fface->front_edges[small_angle_idx0];
    main_fedges[1] = fface->front_edges[small_angle_idx1];
    std::array<Edge*, 2> main_edges;
    main_edges[0] = main_fedges[0]->x;
    main_edges[1] = main_fedges[1]->x;
    auto main_vert = main_edges[0]->verts[0] == main_edges[1]->verts[0] ?
        main_edges[0]->verts[0] :
        (main_edges[0]->verts[0] == main_edges[1]->verts[1] ?
         main_edges[0]->verts[0] :
         main_edges[0]->verts[1]);

    auto v0 = main_edges[0]->verts[0] == main_vert ? main_edges[0]->verts[1] : main_edges[0]->verts[0];
    auto v1 = main_edges[1]->verts[0] == main_vert ? main_edges[1]->verts[1] : main_edges[1]->verts[0];
    auto v2 = main_vert;

    auto adj_ffaces = main_fedges[0]->adj_ffaces();
    auto fn0 = std::get<0>(adj_ffaces) == fface ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn0 = fn0->x->find_vert_not(main_fedges[0]->x);
    adj_ffaces = main_fedges[1]->adj_ffaces();
    auto fn1 = std::get<0>(adj_ffaces) == fface ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn1 = fn1->x->find_vert_not(main_fedges[1]->x);

    vec3 v0_pos = v0->pos();
    vec3 v1_pos = v1->pos();
    vec3 v2_pos = v2->pos();
    vec3 vn0_pos = vn0->pos();
    vec3 vn1_pos = vn1->pos();

    vec3 n_m = fface->normal;
    vec3 n_n0 = fn0->normal;
    vec3 n_n1 = fn1->normal;

    vec3 e_mn0 = (n_m + n_n0).normalize();
    vec3 e_mn1 = (n_m + n_n1).normalize();

    vec3 np_mn0 = spt::cross(v2_pos - v0_pos, e_mn0);
    vec3 np_mn1 = spt::cross(v2_pos - v1_pos, e_mn1);

    vec3 e = spt::cross(np_mn0, np_mn1).normalize();
    if (spt::dot(e, n_m) < static_cast<real_t>(0.0)) e *= static_cast<real_t>(-1.0);

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
    vec3 new_pos = v2_pos + e * magn_d;

    if (edge_intersect_front(v0, new_pos) ||
        edge_intersect_front(v1, new_pos) ||
        edge_intersect_front(v2, new_pos) ||
        does_front_intersect_sphere(new_pos, C_MIN_DIS * magn_d) ||
        any_edge_intersect_face(v0, v1, new_pos) ||
        any_edge_intersect_face(v0, v2, new_pos) ||
        any_edge_intersect_face(v1, v2, new_pos)) {
        // NOTE: do i really need it?

        new_pos = new_pos = v2_pos + e * av_magn;
        if (edge_intersect_front(v0, new_pos) ||
            edge_intersect_front(v1, new_pos) ||
            edge_intersect_front(v2, new_pos) ||
            does_front_intersect_sphere(new_pos, C_MIN_DIS * av_magn) ||
            any_edge_intersect_face(v0, v1, new_pos) ||
            any_edge_intersect_face(v0, v2, new_pos) ||
            any_edge_intersect_face(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type2";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::try_compute_new_vert_pos_type1(front::Face* fface, vec3& out_pos, std::size_t small_angle_idx) {
    auto main_f_edge = fface->front_edges[small_angle_idx];
    auto main_edge = fface->front_edges[small_angle_idx]->x;
    std::array<vec3, 2> main_edge_poses;
    main_edge_poses[0] = main_edge->verts[0]->pos();
    main_edge_poses[1] = main_edge->verts[1]->pos();

    auto v0 = main_edge->verts[0];
    auto v1 = main_edge->verts[1];
    auto v2 = fface->x->find_vert_not(main_edge);

    auto adj_ffaces = main_f_edge->adj_ffaces();
    auto fn = std::get<0>(adj_ffaces) == fface ? std::get<1>(adj_ffaces) : std::get<0>(adj_ffaces);
    auto vn = fn->x->find_vert_not(main_edge);

    vec3 v0_pos = v0->pos();
    vec3 v1_pos = v1->pos();
    vec3 v2_pos = v2->pos();
    vec3 vn_pos = vn->pos();

    vec3 v2pr = spt::project(v2_pos, v0_pos, v1_pos);
    vec3 vnpr = spt::project(vn_pos, v0_pos, v1_pos);

    vec3 c = (v0_pos + v1_pos + v2pr + vnpr) * static_cast<real_t>(0.25);
    vec3 e = (fface->normal + fn->normal).normalize();
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
    vec3 new_pos = c + e * std::sqrt(magn_d * magn_d - v0_c_dist * v0_c_dist);

    if (edge_intersect_front(v0, new_pos) ||
        edge_intersect_front(v1, new_pos) ||
        edge_intersect_front(v2, new_pos) ||
        does_front_intersect_sphere(new_pos, C_MIN_DIS * magn_d) ||
        any_edge_intersect_face(v0, v1, new_pos) ||
        any_edge_intersect_face(v0, v2, new_pos) ||
        any_edge_intersect_face(v1, v2, new_pos)) {
        // NOTE: do i really need it?

        new_pos = c + e * std::sqrt(av_magn * av_magn - v0_c_dist * v0_c_dist);
        if (edge_intersect_front(v0, new_pos) ||
            edge_intersect_front(v1, new_pos) ||
            edge_intersect_front(v2, new_pos) ||
            does_front_intersect_sphere(new_pos, C_MIN_DIS * av_magn) ||
            any_edge_intersect_face(v0, v1, new_pos) ||
            any_edge_intersect_face(v0, v2, new_pos) ||
            any_edge_intersect_face(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type1";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::try_compute_new_vert_pos_type0(front::Face* fface, vec3& out_pos) {
    real_t av_magn = ONE_3 * (fface->x->edges[0]->magnitude()
                              + fface->x->edges[1]->magnitude()
                              + fface->x->edges[2]->magnitude());
    real_t raw_deform = C_D * (m_prefLen - av_magn);
    real_t deform = raw_deform < av_magn * C_MAXD ? raw_deform : av_magn * C_MAXD;
    real_t magn_d = av_magn + deform;
    vec3 new_pos = fface->center() + fface->normal * std::sqrt(magn_d * magn_d - ONE_3 * av_magn * av_magn);

    auto v0 = fface->x->edges[0]->verts[0];
    auto v1 = fface->x->edges[0]->verts[1];
    auto v2 = fface->x->find_vert_not(fface->x->edges[0]);

    if (edge_intersect_front(v0, new_pos) ||
        edge_intersect_front(v1, new_pos) ||
        edge_intersect_front(v2, new_pos) ||
        does_front_intersect_sphere(new_pos, C_MIN_DIS * magn_d) ||
        any_edge_intersect_face(v0, v1, new_pos) ||
        any_edge_intersect_face(v0, v2, new_pos) ||
        any_edge_intersect_face(v1, v2, new_pos)) {
        // NOTE: do i really need it?

        new_pos = fface->center() + fface->normal * std::sqrt(av_magn * av_magn - ONE_3 * av_magn * av_magn);
        if (edge_intersect_front(v0, new_pos) ||
            edge_intersect_front(v1, new_pos) ||
            edge_intersect_front(v2, new_pos) ||
            does_front_intersect_sphere(new_pos, C_MIN_DIS * av_magn) ||
            any_edge_intersect_face(v0, v1, new_pos) ||
            any_edge_intersect_face(v0, v2, new_pos) ||
            any_edge_intersect_face(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type0";

    out_pos = new_pos;
    return true;
}


bool Polyhedron::try_compute_new_vert_pos(front::Face* fface, vec3& out_pos) {
    std::array<real_t, 3> angs {
        fface->front_edges[0]->angle(),
        fface->front_edges[1]->angle(),
        fface->front_edges[2]->angle()
    };
    std::array<std::size_t, 3> idces;
    std::size_t n_small_angs = 0;
    if (angs[0] < degToRad(140)) idces[n_small_angs++] = 0;
    if (angs[1] < degToRad(140)) idces[n_small_angs++] = 1;
    if (angs[2] < degToRad(140)) idces[n_small_angs++] = 2;

    switch (n_small_angs) {
    case 0: return try_compute_new_vert_pos_type0(fface, out_pos);
    case 1: return try_compute_new_vert_pos_type1(fface, out_pos, idces[0]);
    case 2: return try_compute_new_vert_pos_type2(fface, out_pos, idces[0], idces[1]);
    case 3: return try_compute_new_vert_pos_type3(fface, out_pos);
    }

    return true;
}


real_t Polyhedron::sqr_face_4areas(const front::Face* fface) const {
    vec3 edge0_vec = fface->x->edges[0]->verts[1]->pos() - fface->x->edges[0]->verts[0]->pos();
    vec3 edge1_vec = fface->x->edges[1]->verts[1]->pos() - fface->x->edges[1]->verts[0]->pos();
    return spt::cross(edge0_vec, edge1_vec).sqr_magnitude();
}


front::Face* Polyhedron::choose_fface_for_exhaustion_with_new_vert(front::Edge* fedge) {
    auto adj_ffaces = fedge->adj_ffaces();

    return sqr_face_4areas(std::get<0>(adj_ffaces)) < sqr_face_4areas(std::get<1>(adj_ffaces)) ?
//    return std::get<0>(adj_ffaces)->quality() < std::get<1>(adj_ffaces)->quality() ?
        std::get<0>(adj_ffaces) : std::get<1>(adj_ffaces);
}


void Polyhedron::exhaust_with_new_vert(front::Face* fface, const vec3& vert_pos) {
    auto new_fvert = add_to_front(new Vert(vert_pos));

    std::array<front::Edge*, 6> new_tetr_fedges;
    new_tetr_fedges[0] = fface->front_edges[0];
    auto far_fvert = fface->find_front_vert_not(new_tetr_fedges[0]);
    new_tetr_fedges[1] = fface->find_front_edge(new_tetr_fedges[0]->front_verts[1], far_fvert);
    new_tetr_fedges[2] = fface->find_front_edge(new_tetr_fedges[0]->front_verts[0], far_fvert);
    new_tetr_fedges[3] = add_to_front(new front::Edge(this, new_tetr_fedges[0]->front_verts[0], new_fvert));
    new_tetr_fedges[4] = add_to_front(new front::Edge(this, new_tetr_fedges[0]->front_verts[1], new_fvert));
    new_tetr_fedges[5] = add_to_front(new front::Edge(this, far_fvert, new_fvert));

    new_tetr_fedges[0]->refresh_angle_data();
    new_tetr_fedges[1]->refresh_angle_data();
    new_tetr_fedges[2]->refresh_angle_data();

    new_tetr_fedges[0]->remove_adj_fface(fface);
    new_tetr_fedges[1]->remove_adj_fface(fface);
    new_tetr_fedges[2]->remove_adj_fface(fface);

    auto new_fface = add_to_front(new front::Face(this, new_tetr_fedges[0],
                                                        new_tetr_fedges[3],
                                                        new_tetr_fedges[4]));
    new_fface->normal = compute_normal_in_tetr(new_fface, far_fvert->x->pos());

    new_fface = add_to_front(new front::Face(this, new_tetr_fedges[2],
                                                   new_tetr_fedges[3],
                                                   new_tetr_fedges[5]));
    new_fface->normal = compute_normal_in_tetr(new_fface, new_tetr_fedges[0]->front_verts[1]->x->pos());

    new_fface = add_to_front(new front::Face(this, new_tetr_fedges[1],
                                                   new_tetr_fedges[4],
                                                   new_tetr_fedges[5]));
    new_fface->normal = compute_normal_in_tetr(new_fface, new_tetr_fedges[0]->front_verts[0]->x->pos());
    
    m_inner_tetrs.push_back(new Tetr(new_tetr_fedges[0]->front_verts[0]->x,
                                     new_tetr_fedges[0]->front_verts[1]->x,
                                     new_tetr_fedges[5]->front_verts[0]->x,
                                     new_tetr_fedges[5]->front_verts[1]->x));

    remove_from_front(fface);
    delete fface;

#ifdef DEBUG
    debug();
    std::cout << std::endl << "exhaust_with_new_vert";
#endif // DEBUG
}


bool Polyhedron::try_exhaust_without_new_vert(front::Edge* fedge) {
    // TODO: improve that checks
    if (will_parallel_faces(fedge) || // TODO: replace later with check if front is close to potential element
        will_any_edge_intersect_faces(fedge) ||
        will_any_vert_inside_tetr(fedge))
        return false;

    exhaust_without_new_vert(fedge);
    return true;
}


bool Polyhedron::try_exhaust_with_new_vert(front::Edge* frontEdge) {
    // TODO: replace later with check if front is close to potential element
    if (will_parallel_faces(frontEdge))
        return false;

    auto exhaust_fface = choose_fface_for_exhaustion_with_new_vert(frontEdge);
    vec3 new_vert_pos;
    if (!try_compute_new_vert_pos(exhaust_fface, new_vert_pos))
        return false;

    exhaust_with_new_vert(exhaust_fface, new_vert_pos);
    return true;
}


bool Polyhedron::front_exhausted() {
    if (m_front_faces.size() == 0 &&
        m_front_edges.size() == 0)
        return true;

    if ((m_front_faces.size() == 0 && m_front_edges.size() > 0) ||
        (m_front_faces.size() > 0 && m_front_edges.size() == 0))
        throw std::logic_error("Error in Polyhedron::front_exhausted. Front wasn't correctly exhausted.");

    return false;
}


void Polyhedron::process_angles() {
#ifdef DEBUG
    int debug_i = 0;
#endif // DEBUG
    real_t max_compl = std::numeric_limits<real_t>::max();
    for (front::Edge* cur_fedge = current_front_edge(max_compl);; cur_fedge = current_front_edge(max_compl)) {
        if (!cur_fedge)
            throw std::logic_error("pmg::Polyhedron::current_front_edge returned nullptr");

        if (exhaust_without_new_vert_priority_predicate(cur_fedge)) {
            if (!try_exhaust_without_new_vert(cur_fedge)) {
                max_compl = cur_fedge->complexity();
                continue;
            }
        } else if (exhaust_with_new_vert_priority_predicate(cur_fedge)) {
            if (!try_exhaust_with_new_vert(cur_fedge)) {
                max_compl = cur_fedge->complexity();
                continue;
            }
        } else {
            front::Face* exhaust_from_fface = nullptr;
            vec3* new_vert_pos = nullptr;
            switch (exhaustion_type_quality_priority(cur_fedge, exhaust_from_fface, new_vert_pos)) {
            case exhaust_type::without_new_vert:
                exhaust_without_new_vert(cur_fedge);
                break;

            case exhaust_type::with_new_vert:
                if (new_vert_pos) {
                    exhaust_with_new_vert(exhaust_from_fface, *new_vert_pos);
                    delete new_vert_pos;
                } else {
                    if (!try_exhaust_with_new_vert(cur_fedge)) {
                        max_compl = cur_fedge->complexity();
                        continue;
                    }
                }
                break;

            case exhaust_type::dont_exhaust:
                max_compl = cur_fedge->complexity();
#ifdef DEBUG
                std::cout << std::endl << "dont_exhaust";
#endif // DEBUG
                continue;
            }
        }
        max_compl = std::numeric_limits<real_t>::max();
        
#ifdef DEBUG
        std::cout << std::endl << debug_i;
        if (debug_i++ == 159)
            m_polyhset->output(filetype::wavefront_obj, "debug.obj");
        debug();
#endif // DEBUG

        if (front_exhausted())
            return;
    }
}

#ifdef DEBUG
void Polyhedron::debug() {
    real_t accum = static_cast<real_t>(0);
    for (auto& fface : m_front_faces) {
        accum += fface->front_edges[0]->angle();
        accum += fface->front_edges[1]->angle();
        accum += fface->front_edges[2]->angle();
    }

    real_t sqrlen = static_cast<real_t>(0);
    for (auto& fface : m_front_faces) {
        accum += fface->front_edges[0]->x->sqr_magnitude();
        accum += fface->front_edges[1]->x->sqr_magnitude();
        accum += fface->front_edges[2]->x->sqr_magnitude();
    }
    
    std::cout << "\n{ " << accum + sqrlen << " }";

    for (auto& fedge : m_front_edges) {
        std::size_t count = 0;
        for (auto& fface : m_front_faces)
            if (fface->contains(fedge))
                count++;

        if (count != 2)
            throw std::exception();
    }

    for (auto& fface : m_front_faces)
        for (auto& fedge : fface->front_edges)
            if (!fface->x->contains(fedge->x))
                throw std::exception();

    for (auto& fface : m_front_faces)
        for (auto& edge : fface->x->edges)
            if (!fface->find_front_edge(edge))
                throw std::exception();
}
#endif // DEBUG


void Polyhedron::tetrahedralize(real_t preferredLen, genparams::Volume gen_params) {
    m_prefLen = preferredLen;
    initialize_front();
    compute_front_normals();
    process_angles();
//    if (global_intersection())
//        throw std::logic_error("Intersection error.\npmg::Polyhedron::global_intersection returned true.");

    smooth_mesh(gen_params.nSmoothIters);
}


bool Polyhedron::global_intersection() {
    for (auto& edge : m_inner_edges)
        if (edge_intersect_any_face(edge))
            return true;

    return false;
}


void Polyhedron::smooth_mesh(std::size_t nIters) {
    for (std::size_t i = 0; i < nIters; i++) {
        for (auto& vert : m_inner_verts) {
            vec3 shift;
            std::size_t n_delta_shifts = 0;
            for (auto& edge : m_inner_edges) {
                if (vert == edge->verts[0]) {
                    shift += edge->verts[1]->pos() - vert->pos();
                    n_delta_shifts++;
                } else if (vert == edge->verts[1]) {
                    shift += edge->verts[0]->pos() - vert->pos();
                    n_delta_shifts++;
                }
            }
            shift /= n_delta_shifts;
            vert->pos() = vert->pos() + shift;
        }
    }
}


pair_rr Polyhedron::analyze_mesh_quality(std::vector<Tetr*>::iterator* out_min_quality_tetr) {
    std::size_t n_tetrs = 0;
    real_t av_q = static_cast<real_t>(0);
    real_t min_q = static_cast<real_t>(1);
    auto min_q_tetr = m_inner_tetrs.begin();
    for (auto iter = m_inner_tetrs.begin(); iter != m_inner_tetrs.end(); iter++) {
        real_t q = (*iter)->quality();
        av_q += q;
        n_tetrs++;
        if (q < min_q) {
            min_q = q;
            if (out_min_quality_tetr)
                min_q_tetr = iter;
        }
    }
    av_q /= n_tetrs;
    if (out_min_quality_tetr)
        *out_min_quality_tetr = min_q_tetr;

    return m_meshQuality = { min_q, av_q };
}


pair_rr Polyhedron::analyze_mesh_abs_grad() {
    for (auto& tetr : m_inner_tetrs) {
        real_t volume = tetr->volume();
        for (auto& vert : tetr->verts) {
            vert->min_adj_tetr_vol = std::min(volume, vert->min_adj_tetr_vol);
            vert->max_adj_tetr_vol = std::max(volume, vert->max_adj_tetr_vol);
        }
    }

    real_t min_abs_grad = static_cast<real_t>(1.0);
    real_t av_abs_grad = static_cast<real_t>(0.0);
    for (auto& vert : m_inner_verts) {
        real_t cur_abs_grad = vert->min_adj_tetr_vol / vert->max_adj_tetr_vol;
        min_abs_grad = std::min(cur_abs_grad, min_abs_grad);
        av_abs_grad += cur_abs_grad;
    }
    av_abs_grad /= m_inner_verts.size();

    return m_meshAbsGrad = { min_abs_grad, av_abs_grad };
}




Polyhedron::Polyhedron(PolyhedralSet* polyhset) {
    m_polyhset = polyhset;
}


Polyhedron::~Polyhedron() {
    for (auto& tetr : m_inner_tetrs)
        delete tetr;
    for (auto& face : m_inner_faces)
        delete face;
    for (auto& edge : m_inner_edges)
        delete edge;
    for (auto& vert : m_inner_verts)
        delete vert;
}




real_t Polyhedron::preferred_length() const {
    return m_prefLen;
}


void Polyhedron::add_to_shell(const shell::Face* shellFace) {
    m_shell_faces.push_back(const_cast<shell::Face*>(shellFace));
}


void Polyhedron::add_to_shell(const shell::Edge* shellEdge) {
    m_shell_edges.push_back(const_cast<shell::Edge*>(shellEdge));
}


void Polyhedron::add_to_shell(const shell::Vert* shellVert) {
    m_shell_verts.push_back(const_cast<shell::Vert*>(shellVert));
}


bool Polyhedron::shell_contains(const shell::Face* shellFace) const {
    return std::find(m_shell_faces.begin(), m_shell_faces.end(), shellFace) != m_shell_faces.end();
}


bool Polyhedron::shell_contains(const shell::Edge* shellEdge) const {
    return std::find(m_shell_edges.begin(), m_shell_edges.end(), shellEdge) != m_shell_edges.end();
}


bool Polyhedron::shell_contains(const shell::Vert* shellVert) const {
    return std::find(m_shell_verts.begin(), m_shell_verts.end(), shellVert) != m_shell_verts.end();
}


const std::vector<Tetr*>& Polyhedron::inner_tetrs() const {
    return m_inner_tetrs;
}


const std::vector<Face*>& Polyhedron::inner_faces() const {
    return m_inner_faces;
}


const std::vector<Vert*>& Polyhedron::inner_verts() const {
    return m_inner_verts;
}


const std::vector<front::Face*>& Polyhedron::front_faces() const {
    return m_front_faces;
}


const std::vector<front::Edge*>& Polyhedron::front_edges() const {
    return m_front_edges;
}
