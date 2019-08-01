// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

// With help of David Eberly's and Dan Sunday's works.

#pragma once
#include <cmath>
#include <algorithm>
#include "mat.h"

// TODO: try pass by classes-containers instead of reference and make benchmark
// TODO: try vectorization (SIMD) and make benchmark
// TODO: try refactor functions
// https://en.cppreference.com/w/cpp/types/numeric_limits/std::numeric_limits<T>::epsilon()
namespace spt {

template <typename T>
bool weak_between(T boundary0, T boundary1, T value) {
    auto l_eps = std::numeric_limits<T>::epsilon() * std::abs(value);
    return
        (value > std::min(boundary0, boundary1) - l_eps) &&
        (value < std::max(boundary0, boundary1) + l_eps);
}

template <typename T>
bool weak_in_cuboid(T corner0, T corner1, T point) {
    return
        weak_between(corner0.x[0], corner1.x[0], point.x[0]) &&
        weak_between(corner0.x[1], corner1.x[1], point.x[1]) &&
        weak_between(corner0.x[2], corner1.x[2], point.x[2]);
}

template <typename Real>
mat<3, Real> dot(const mat<3, Real>& mat0, const mat<3, Real>& mat1) {
    std::array mat1_cols{
        vec<3, Real>{ mat1[0][0], mat1[1][0], mat1[2][0] },
        vec<3, Real>{ mat1[0][1], mat1[1][1], mat1[2][1] },
        vec<3, Real>{ mat1[0][2], mat1[1][2], mat1[2][2] } };
    return {
        { dot(mat0[0], mat1_cols[0]), dot(mat0[0], mat1_cols[1]), dot(mat0[0], mat1_cols[2]) },
        { dot(mat0[1], mat1_cols[0]), dot(mat0[1], mat1_cols[1]), dot(mat0[1], mat1_cols[2]) },
        { dot(mat0[2], mat1_cols[0]), dot(mat0[2], mat1_cols[1]), dot(mat0[2], mat1_cols[2]) } };
}

template <typename Real>
vec<3, Real> dot(const mat<3, Real>& matr, const vec<3, Real>& vect) {
    return { dot(matr[0], vect), dot(matr[1], vect), dot(matr[2], vect) };
}

template <typename Real>
Real dot(const vec<3, Real>& vec0, const vec<3, Real>& vec1) {
    return vec0[0] * vec1[0] + vec0[1] * vec1[1] + vec0[2] * vec1[2];
}

template <typename Real>
Real dot(const vec<2, Real>& vec0, const vec<2, Real>& vec1) {
    return vec0[0] * vec1[0] + vec0[1] * vec1[1];
}

template <typename Real>
vec<3, Real> cross(const vec<3, Real>& vec0, const vec<3, Real>& vec1) {
    return {
        vec0.x[1] * vec1.x[2] - vec0.x[2] * vec1.x[1],
        vec0.x[2] * vec1.x[0] - vec0.x[0] * vec1.x[2],
        vec0.x[0] * vec1.x[1] - vec0.x[1] * vec1.x[0] };
}

template <typename Real>
Real cross(const vec<2, Real>& vec0, const vec<2, Real>& vec1) {
    return vec0[0] * vec1[1] - vec0[1] * vec1[0];
}

template <typename Real>
Real mixed(const vec<3, Real>& vec0, const vec<3, Real>& vec1, const vec<3, Real>& vec2) {
    return dot(cross(vec0, vec1), vec2);
}

template <typename Real>
Real cos(const vec<3, Real>& vec0, const vec<3, Real>& vec1) {
    return dot(vec0, vec1) / std::sqrt(vec0.sqr_magnitude() * vec1.sqr_magnitude());
}

template <typename Real>
Real cos(const vec<2, Real>& vec0, const vec<2, Real>& vec1) {
    return dot(vec0, vec1) / std::sqrt(vec0.sqr_magnitude() * vec1.sqr_magnitude());
}

template <typename Real>
vec<3, Real> project(
    const vec<3, Real>& point,
    const vec<3, Real>& line_p0, const vec<3, Real>& line_p1) {

    return line_p0 + (point - line_p0).project(line_p1 - line_p0);
}

template <typename Real>
bool project(
    vec<3, Real>& out,
    const vec<3, Real>& point,
    const vec<3, Real>& segm_p0, const vec<3, Real>& segm_p1) {

    auto res = segm_p0 + (point - segm_p0).project(segm_p1 - segm_p0);

    if (in_cuboid(segm_p0, segm_p1, res)) {
        out = res;
        return true;
    }

    return false;
}

template <typename Real>
vec<3, Real> project(
    const vec<3, Real>& point,
    const vec<3, Real>& plane_p0, const vec<3, Real>& plane_p1, const vec<3, Real>& plane_p2) {

    return plane_p0 + (point - plane_p0).project(plane_p1 - plane_p0, plane_p2 - plane_p0);
}

template <typename Real>
bool does_ray_intersect_plane(
    const vec<3, Real>& dir,
    const vec<3, Real>& pl_p0, const vec<3, Real>& pl_p1, const vec<3, Real>& pl_p2) {

    std::array edges{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    return std::abs(det) > std::numeric_limits<Real>::epsilon();
}

// TODO: return std::optional<vec<3>>
template <typename Real>
bool ray_intersect_plane(
    vec<3, Real>& out_intersect_point,
    const vec<3, Real>& origin, const vec<3, Real>& dir,
    const vec<3, Real>& pl_p0, const vec<3, Real>& pl_p1, const vec<3, Real>& pl_p2) {

    std::array edges{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto tvec = origin - pl_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersect_point = origin + dir * t;
    return t > static_cast<Real>(0);
}

template <typename Real>
bool does_ray_intersect_triangle(
    const vec<3, Real>& origin, const vec<3, Real>& dir,
    const vec<3, Real>& tr_p0, const vec<3, Real>& tr_p1, const vec<3, Real>& tr_p2) {

    std::array edges{ tr_p1 - tr_p0, tr_p2 - tr_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto inv_det = static_cast<Real>(1) / det;

    auto tvec = origin - tr_p0;
    auto u = spt::dot(tvec, pvec) * inv_det;
    if (u < static_cast<Real>(0) || u > static_cast<Real>(1))
        return false;

    auto qvec = spt::cross(tvec, edges[0]);
    auto v = spt::dot(dir, qvec) * inv_det;
    if (v < static_cast<Real>(0) || u + v > static_cast<Real>(1))
        return false;

    auto t = spt::dot(edges[1], qvec) * inv_det;
    return t >= static_cast<Real>(0);
}

// TODO: return std::optional<vec<3>>
template <typename Real>
bool line_intersect_plane(
    vec<3, Real>& out_intersect_point,
    const vec<3, Real>& line_point, const vec<3, Real>& line_dir,
    const vec<3, Real>& plane_p0, const vec<3, Real>& plane_p1, const vec<3, Real>& plane_p2) {

    std::array edges{ plane_p1 - plane_p0, plane_p2 - plane_p0 };

    auto pvec = spt::cross(line_dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto tvec = line_point - plane_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersect_point = line_point + line_dir * t;
    return true;
}

template <typename Real>
vec<3, Real> line_intersect_plane(
    const vec<3, Real>& line_point, const vec<3, Real>& line_dir,
    const vec<3, Real>& plane_p0, const vec<3, Real>& plane_p1, const vec<3, Real>& plane_p2) {

    std::array edges{ plane_p1 - plane_p0, plane_p2 - plane_p0 };

    auto pvec = spt::cross(line_dir, edges[1]);
    auto tvec = line_point - plane_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / spt::dot(edges[0], pvec);
    return line_point + line_dir * t;
}

template <typename Real>
bool does_segment_intersect_triangle(
    const vec<3, Real>& segm_p0, const vec<3, Real>& segm_p1,
    const vec<3, Real>& tr_p0, const vec<3, Real>& tr_p1, const vec<3, Real>& tr_p2) {

    const auto dir = segm_p1 - segm_p0;

    std::array edges{ tr_p1 - tr_p0, tr_p2 - tr_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto inv_det = static_cast<Real>(1) / det;

    auto tvec = segm_p0 - tr_p0;
    auto u = spt::dot(tvec, pvec) * inv_det;
    if (u < static_cast<Real>(0) || u > static_cast<Real>(1))
        return false;

    auto qvec = spt::cross(tvec, edges[0]);
    auto v = spt::dot(dir, qvec) * inv_det;
    if (v < static_cast<Real>(0) || u + v > static_cast<Real>(1))
        return false;

    auto t = spt::dot(edges[1], qvec) * inv_det;
    return t <= static_cast<Real>(1) && t >= static_cast<Real>(0);
}

// TODO: return std::optional<vec<3>>
template <typename Real>
bool segment_intersect_plane(
    vec<3, Real>& out_intersect_point,
    const vec<3, Real>& p0, const vec<3, Real>& p1,
    const vec<3, Real>& pl_p0, const vec<3, Real>& pl_p1, const vec<3, Real>& pl_p2) {

    const auto dir = p1 - p0;

    std::array edges{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (std::abs(det) <= std::numeric_limits<Real>::epsilon())
        return false;

    auto tvec = p0 - pl_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersect_point = p0 + dir * t;
    return t <= static_cast<Real>(1) && t >= static_cast<Real>(0);
}

template <typename Real>
bool does_triangle_intersect_sphere(
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2,
    const vec<3, Real>& center, Real radius) {

    auto proj = project(center, trngl_p0, trngl_p1, trngl_p2);
    if ((proj - center).sqr_magnitude() > radius * radius)
        return false;

    if (is_point_on_triangle(proj, trngl_p0, trngl_p1, trngl_p2, max_sqrs_sum(trngl_p0, trngl_p1, trngl_p2)))
        return true;

    auto closest = closest_triangle_point_to_point_on_plane(proj, trngl_p0, trngl_p1, trngl_p2);
    return (closest - center).sqr_magnitude() <= radius * radius;
}

template <typename Real>
Real sqrs_sum(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    std::array sqrs{
        (trngl_p0 - point).sqr_magnitude(),
        (trngl_p1 - point).sqr_magnitude(),
        (trngl_p2 - point).sqr_magnitude()
    };
    return sqrs[0] + sqrs[1] + sqrs[2];
}

template <typename Real>
Real max_sqrs_sum(
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    std::array sqrs{
        (trngl_p1 - trngl_p0).sqr_magnitude(),
        (trngl_p2 - trngl_p1).sqr_magnitude(),
        (trngl_p0 - trngl_p2).sqr_magnitude()
    };

    std::array<std::size_t, 2> max_inds;
    if (sqrs[0] < sqrs[1]) {
        max_inds[0] = 1;
        if (sqrs[0] < sqrs[2])
            max_inds[1] = 2;
        else
            max_inds[1] = 0;
    } else {
        max_inds[0] = 0;
        if (sqrs[1] < sqrs[2])
            max_inds[1] = 2;
        else
            max_inds[1] = 1;
    }

    return sqrs[max_inds[0]] + sqrs[max_inds[1]];
}

template <typename Real>
bool is_point_on_triangle(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    auto s0 = spt::cross(trngl_p0 - point, trngl_p1 - point).magnitude();
    auto s1 = spt::cross(trngl_p0 - point, trngl_p2 - point).magnitude();
    auto s2 = spt::cross(trngl_p1 - point, trngl_p2 - point).magnitude();
    auto s = spt::cross(trngl_p0 - trngl_p2, trngl_p1 - trngl_p2).magnitude();

    auto expr = s - s0 - s1 - s2;
    return std::abs(expr) <= std::numeric_limits<Real>::epsilon();
}

template <typename Real>
bool is_point_on_triangle(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2,
    Real max_sqrs_sum) {

    if (sqrs_sum(point, trngl_p0, trngl_p1, trngl_p2) > max_sqrs_sum)
        return false;

    return is_point_on_triangle(point, trngl_p0, trngl_p1, trngl_p2);
}

template <typename Real>
bool is_point_in_tetrahedron(
    const vec<3, Real>& point,
    const vec<3, Real>& tetr_p0, const vec<3, Real>& tetr_p1,
    const vec<3, Real>& tetr_p2, const vec<3, Real>& tetr_p3) {

    auto vert_to_p0 = tetr_p0 - point;
    auto vert_to_p1 = tetr_p1 - point;
    auto vert_to_p2 = tetr_p2 - point;
    auto vert_to_p3 = tetr_p3 - point;

    std::array abs_mixed_prods{
        std::abs(spt::mixed(vert_to_p0, vert_to_p2, vert_to_p3)),
        std::abs(spt::mixed(vert_to_p0, vert_to_p1, vert_to_p2)),
        std::abs(spt::mixed(vert_to_p0, vert_to_p1, vert_to_p3)),
        std::abs(spt::mixed(vert_to_p1, vert_to_p2, vert_to_p3)),
        std::abs(spt::mixed(tetr_p1 - tetr_p0, tetr_p2 - tetr_p0, tetr_p3 - tetr_p0))
    };

    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4];
}

template <typename Real>
vec<3, Real> closest_segment_point_to_point(
    const vec<3, Real>& point,
    const vec<3, Real>& segm_p0, const vec<3, Real>& segm_p1) {

    auto proj = project(point, segm_p0, segm_p1);

    if (weak_in_cuboid(segm_p0, segm_p1, proj)) {
        return proj;
    } else if (std::array sqr_magns{ (segm_p0 - point).sqr_magnitude(), (segm_p1 - point).sqr_magnitude() };
               sqr_magns[0] < sqr_magns[1]) {
        return segm_p0;
    } else {
        return segm_p1;
    }
}

template <typename Real>
vec<3, Real> closest_triangle_point_to_point_on_plane(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    std::array closest_points{
        closest_segment_point_to_point(point, trngl_p0, trngl_p1),
        closest_segment_point_to_point(point, trngl_p1, trngl_p2),
        closest_segment_point_to_point(point, trngl_p2, trngl_p0)
    };

    std::array sqrs{
        (closest_points[0] - point).sqr_magnitude(),
        (closest_points[1] - point).sqr_magnitude(),
        (closest_points[2] - point).sqr_magnitude()
    };

    std::size_t min_i;
    if (sqrs[0] < sqrs[1]) {
        if (sqrs[0] < sqrs[2])
            min_i = 0;
        else
            min_i = 2;
    } else {
        if (sqrs[1] < sqrs[2])
            min_i = 1;
        else
            min_i = 2;
    }

    return closest_points[min_i];
}

template <typename Real>
Real distance_point_to_line(
    const vec<3, Real>& point,
    const vec<3, Real>& line_p0, const vec<3, Real>& line_p1) {

    return (project(point, line_p0, line_p1) - point).magnitude();
}

template <typename Real>
Real distance_point_to_segment(
    const vec<3, Real>& point,
    const vec<3, Real>& segm_p0, const vec<3, Real>& segm_p1) {

    auto proj = project(point, segm_p0, segm_p1);

    if (weak_in_cuboid(segm_p0, segm_p1, proj)) {
        return (proj - point).magnitude();
    } else if (std::array sqr_magns{ (segm_p0 - point).sqr_magnitude(), (segm_p1 - point).sqr_magnitude() };
               sqr_magns[0] < sqr_magns[1]) {
        return std::sqrt(sqr_magns[0]);
    } else {
        return std::sqrt(sqr_magns[1]);
    }
}

template <typename Real>
Real distance_point_to_triangle_on_plane(
    const vec<3, Real>& point,
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    std::array closest_points{
        closest_segment_point_to_point(point, trngl_p0, trngl_p1),
        closest_segment_point_to_point(point, trngl_p1, trngl_p2),
        closest_segment_point_to_point(point, trngl_p2, trngl_p0)
    };

    std::array sqrs{
        (closest_points[0] - point).sqr_magnitude(),
        (closest_points[1] - point).sqr_magnitude(),
        (closest_points[2] - point).sqr_magnitude()
    };

    return std::min({ sqrs[0], sqrs[1], sqrs[2] });
}

template <typename Real>
std::pair<vec<3, Real>, vec<3, Real>> lines_closest_points(
    const vec<3, Real>& line0_p0, const vec<3, Real>& line0_p1,
    const vec<3, Real>& line1_p0, const vec<3, Real>& line1_p1) {

    auto u = line0_p1 - line0_p0;
    auto v = line1_p1 - line1_p0;
    auto w = line0_p0 - line1_p0;
    auto a = spt::dot(u, u); // always >= 0
    auto b = spt::dot(u, v);
    auto c = spt::dot(v, v); // always >= 0
    auto d = spt::dot(u, w);
    auto e = spt::dot(v, w);
    auto D = a * c - b * b; // always >= 0

    Real sc, sN, sD = D; // sc = sN / sD, default sD = D >= 0
    Real tc, tN, tD = D; // tc = tN / tD, default tD = D >= 0
    // compute the line parameters of the two closest points
    if (D <= std::numeric_limits<Real>::epsilon()) { // the lines are almost parallel
        sN = static_cast<Real>(0); // force using point P0 on segment S1
        sD = static_cast<Real>(1); // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else { // get the closest points on the infinite lines
        sN = b * e - c * d;
        tN = a * e - b * d;
    }

    // finally do the division to get sc and tc
    sc = std::abs(sN) <= std::numeric_limits<Real>::epsilon() ? static_cast<Real>(0) : sN / sD;
    tc = std::abs(tN) <= std::numeric_limits<Real>::epsilon() ? static_cast<Real>(0) : tN / tD;

    return {
        (static_cast<Real>(1) - sc) * line0_p0 + sc * line0_p1,
        (static_cast<Real>(1) - tc) * line1_p0 + tc * line1_p1
    };
}

template <typename Real>
vec<3, Real> lines_closest_point(
    const vec<3, Real>& line0_p0, const vec<3, Real>& line0_p1,
    const vec<3, Real>& line1_p0, const vec<3, Real>& line1_p1) {

    auto closest = lines_closest_points(line0_p0, line0_p1, line1_p0, line1_p1);
    return (std::get<0>(closest) + std::get<1>(closest)) * static_cast<Real>(0.5);
}

template <typename Real>
Real lines_distance(
    const vec<3, Real>& line0_p0, const vec<3, Real>& line0_p1,
    const vec<3, Real>& line1_p0, const vec<3, Real>& line1_p1) {

    auto closest = lines_closest_points(line0_p0, line0_p1, line1_p0, line1_p1);
    auto diff = std::get<0>(closest) - std::get<1>(closest);
    return diff.magnitude();
}

// TODO: pair -> array
template <typename Real>
std::pair<vec<3, Real>, vec<3, Real>> segments_closest_points(
    const vec<3, Real>& segm0_p0, const vec<3, Real>& segm0_p1,
    const vec<3, Real>& segm1_p0, const vec<3, Real>& segm1_p1) {

    auto u = segm0_p1 - segm0_p0;
    auto v = segm1_p1 - segm1_p0;
    auto w = segm0_p0 - segm1_p0;
    auto a = dot(u, u); // always >= 0
    auto b = dot(u, v);
    auto c = dot(v, v); // always >= 0
    auto d = dot(u, w);
    auto e = dot(v, w);
    auto D = a * c - b * b; // always >= 0
    Real sc, sN, sD = D; // sc = sN / sD, default sD = D >= 0
    Real tc, tN, tD = D; // tc = tN / tD, default tD = D >= 0

    // compute the line parameters of the two closest points
    if (D <= std::numeric_limits<Real>::epsilon()) { // the lines are almost parallel
        sN = static_cast<Real>(0); // force using point P0 on segment S1
        sD = static_cast<Real>(1); // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    } else { // get the closest points on the infinite lines
        sN = b * e - c * d;
        tN = a * e - b * d;
        if (sN < static_cast<Real>(0)) { // sc < 0 => the s=0 edge is visible
            sN = static_cast<Real>(0);
            tN = e;
            tD = c;
        } else if (sN > sD) {  // sc > 1  => the s=1 edge is visible
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < static_cast<Real>(0)) { // tc < 0 => the t=0 edge is visible
        tN = static_cast<Real>(0);
        // recompute sc for this edge
        if (-d < static_cast<Real>(0))
            sN = static_cast<Real>(0);
        else if (-d > a)
            sN = sD;
        else {
            sN = -d;
            sD = a;
        }
    } else if (tN > tD) { // tc > 1  => the t=1 edge is visible
        tN = tD;
        // recompute sc for this edge
        if (-d + b < static_cast<Real>(0))
            sN = static_cast<Real>(0);
        else if (-d + b > a)
            sN = sD;
        else {
            sN = -d + b;
            sD = a;
        }
    }
    // finally do the division to get sc and tc
    sc = std::abs(sN) <= std::numeric_limits<Real>::epsilon() ? static_cast<Real>(0) : sN / sD;
    tc = std::abs(tN) <= std::numeric_limits<Real>::epsilon() ? static_cast<Real>(0) : tN / tD;

    return {
        (static_cast<Real>(1) - sc) * segm0_p0 + sc * segm0_p1,
        (static_cast<Real>(1) - tc) * segm1_p0 + tc * segm1_p1
    };
}

template <typename Real>
vec<3, Real> segments_closest_point(
    const vec<3, Real>& segm0_p0, const vec<3, Real>& segm0_p1,
    const vec<3, Real>& segm1_p0, const vec<3, Real>& segm1_p1) {

    auto closest = segments_closest_points(segm0_p0, segm0_p1, segm1_p0, segm1_p1);
    return (std::get<0>(closest) + std::get<1>(closest)) * static_cast<Real>(0.5);
}

template <typename Real>
Real segments_distance(
    const vec<3, Real>& segm0_p0, const vec<3, Real>& segm0_p1,
    const vec<3, Real>& segm1_p0, const vec<3, Real>& segm1_p1) {

    auto closest = segments_closest_points(segm0_p0, segm0_p1, segm1_p0, segm1_p1);
    auto diff = std::get<0>(closest) - std::get<1>(closest);
    return diff.magnitude();
}

// CPA - Closest Point of Approach.
template <typename Real>
Real cpa_time(
    const vec<3, Real>& start0, const vec<3, Real>& vel0,
    const vec<3, Real>& start1, const vec<3, Real>& vel1) {

    auto dv = vel0 - vel1;
    auto dv2 = spt::dot(dv, dv);
    if (dv2 <= std::numeric_limits<Real>::epsilon())
        return static_cast<Real>(0);

    auto w0 = start0 - start1;
    return -spt::dot(w0, dv) / dv2;
}

template <typename Real>
Real cpa_distance(
    const vec<3, Real>& start0, const vec<3, Real>& vel0,
    const vec<3, Real>& start1, const vec<3, Real>& vel1) {

    auto time = cpa_time(start0, vel0, start1, vel1);
    auto p0 = start0 + vel0 * time;
    auto p1 = start1 + vel1 * time;
    return (p1 - p0).magnitude();
}

} // namespace spt
