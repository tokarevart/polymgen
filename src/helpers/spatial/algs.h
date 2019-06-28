// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

// With help of GeomAlgorithms.com website.

#pragma once
#include <cmath>
#include <algorithm>
#include "mat.h"

// TODO: experiment with pass by value instead of reference and make benchmark
namespace spt {

template <typename T>
constexpr T epsilon = T(0);
template <>
constexpr auto epsilon<double> = 1e-6;

namespace helpers {
template <typename T>
T det(T a, T b, T c, T d) {
    return a * d - b * c;
}
template <typename T>
bool between(T boundary0, T boundary1, T value) {
    return
        (value > std::min(boundary0, boundary1) - epsilon<T>) &&
        (value < std::max(boundary0, boundary1) + epsilon<T>);
}
template <typename T>
bool in_rectangle(T corner0, T corner1, T point) {
    return
        between(corner0.x[0], corner1.x[0], point.x[0]) &&
        between(corner0.x[1], corner1.x[1], point.x[1]) &&
        between(corner0.x[2], corner1.x[2], point.x[2]);
}
} // namespace helpers

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
    return vec<3, Real>(
        vec0.x[1] * vec1.x[2] - vec0.x[2] * vec1.x[1],
        vec0.x[2] * vec1.x[0] - vec0.x[0] * vec1.x[2],
        vec0.x[0] * vec1.x[1] - vec0.x[1] * vec1.x[0]);
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

    vec<3, Real> res = segm_p0 + (point - segm_p0).project(segm_p1 - segm_p0);

    if (helpers::in_rectangle(segm_p0, segm_p1, res)) {
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

    using real_t = Real;
    std::array<vec<3, Real>, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    return det > epsilon<real_t> || det < -epsilon<real_t>;
}

// TODO: return std::optional<vec<3>>
template <typename Real>
bool ray_intersect_plane(
    vec<3, Real>& out_intersectPoint,
    const vec<3, Real>& origin, const vec<3, Real>& dir,
    const vec<3, Real>& pl_p0, const vec<3, Real>& pl_p1, const vec<3, Real>& pl_p2) {

    using real_t = Real;
    std::array<vec<3, Real>, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (det < epsilon<real_t> && det > -epsilon<real_t>)
        return false;

    auto tvec = origin - pl_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersectPoint = origin + dir * t;
    return t > static_cast<real_t>(0);
}

template <typename Real>
bool does_ray_intersect_triangle(
    const vec<3, Real>& origin, const vec<3, Real>& dir,
    const vec<3, Real>& tr_p0, const vec<3, Real>& tr_p1, const vec<3, Real>& tr_p2) {

    using real_t = Real;
    std::array<vec<3, Real>, 2> edges = { tr_p1 - tr_p0, tr_p2 - tr_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (det < epsilon<real_t> && det > epsilon<real_t>)
        return false;

    auto inv_det = static_cast<real_t>(1) / det;

    auto tvec = origin - tr_p0;
    auto u = spt::dot(tvec, pvec) * inv_det;
    if (u < static_cast<real_t>(0) || u > static_cast<real_t>(1))
        return false;

    auto qvec = spt::cross(tvec, edges[0]);
    auto v = spt::dot(dir, qvec) * inv_det;
    if (v < static_cast<real_t>(0) || u + v > static_cast<real_t>(1))
        return false;

    auto t = spt::dot(edges[1], qvec) * inv_det;
    return t >= static_cast<real_t>(0);
}

// TODO: return std::optional<vec<3>>
template <typename Real>
bool line_intersect_plane(
    vec<3, Real>& out_intersectPoint,
    const vec<3, Real>& line_point, const vec<3, Real>& line_dir,
    const vec<3, Real>& plane_p0, const vec<3, Real>& plane_p1, const vec<3, Real>& plane_p2) {

    using real_t = Real;
    std::array<vec<3, Real>, 2> edges = { plane_p1 - plane_p0, plane_p2 - plane_p0 };

    auto pvec = spt::cross(line_dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (det < epsilon<real_t> && det > -epsilon<real_t>)
        return false;

    auto tvec = line_point - plane_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersectPoint = line_point + line_dir * t;
    return true;
}

template <typename Real>
vec<3, Real> line_intersect_plane(
    const vec<3, Real>& line_point, const vec<3, Real>& line_dir,
    const vec<3, Real>& plane_p0, const vec<3, Real>& plane_p1, const vec<3, Real>& plane_p2) {

    std::array<vec<3, Real>, 2> edges = { plane_p1 - plane_p0, plane_p2 - plane_p0 };

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

    using real_t = Real;
    const auto dir = segm_p1 - segm_p0;

    std::array<vec<3, Real>, 2> edges = { tr_p1 - tr_p0, tr_p2 - tr_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (det < epsilon<real_t> && det > -epsilon<real_t>)
        return false;

    auto inv_det = static_cast<real_t>(1) / det;

    auto tvec = segm_p0 - tr_p0;
    auto u = spt::dot(tvec, pvec) * inv_det;
    if (u < static_cast<real_t>(0) || u > static_cast<real_t>(1))
        return false;

    auto qvec = spt::cross(tvec, edges[0]);
    auto v = spt::dot(dir, qvec) * inv_det;
    if (v < static_cast<real_t>(0) || u + v > static_cast<real_t>(1))
        return false;

    auto t = spt::dot(edges[1], qvec) * inv_det;
    return t <= static_cast<real_t>(1) && t >= static_cast<real_t>(0);
}

// TODO: return std::optional<vec<3>>
template <typename Real>
bool segment_intersect_plane(
    vec<3, Real>& out_intersectPoint,
    const vec<3, Real>& p0, const vec<3, Real>& p1,
    const vec<3, Real>& pl_p0, const vec<3, Real>& pl_p1, const vec<3, Real>& pl_p2) {

    using real_t = Real;
    const auto dir = p1 - p0;

    std::array<vec<3, Real>, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    auto pvec = spt::cross(dir, edges[1]);
    auto det = spt::dot(edges[0], pvec);
    if (det < epsilon<real_t> && det > -epsilon<real_t>)
        return false;

    auto tvec = p0 - pl_p0;
    auto qvec = spt::cross(tvec, edges[0]);

    auto t = spt::dot(edges[1], qvec) / det;
    out_intersectPoint = p0 + dir * t;
    return t <= static_cast<real_t>(1) && t >= static_cast<real_t>(0);
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

    using real_t = Real;
    std::array<real_t, 3> sqrs;
    sqrs[0] = (trngl_p0 - point).sqr_magnitude();
    sqrs[1] = (trngl_p1 - point).sqr_magnitude();
    sqrs[2] = (trngl_p2 - point).sqr_magnitude();
    return sqrs[0] + sqrs[1] + sqrs[2];
}

template <typename Real>
Real max_sqrs_sum(
    const vec<3, Real>& trngl_p0, const vec<3, Real>& trngl_p1, const vec<3, Real>& trngl_p2) {

    using real_t = Real;
    std::array<real_t, 3> sqrs;
    sqrs[0] = (trngl_p1 - trngl_p0).sqr_magnitude();
    sqrs[1] = (trngl_p2 - trngl_p1).sqr_magnitude();
    sqrs[2] = (trngl_p0 - trngl_p2).sqr_magnitude();

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

    using real_t = Real;
    auto s0 = spt::cross(trngl_p0 - point, trngl_p1 - point).magnitude();
    auto s1 = spt::cross(trngl_p0 - point, trngl_p2 - point).magnitude();
    auto s2 = spt::cross(trngl_p1 - point, trngl_p2 - point).magnitude();
    auto s = spt::cross(trngl_p0 - trngl_p2, trngl_p1 - trngl_p2).magnitude();

    auto expr = s - s0 - s1 - s2;
    return expr > -epsilon<real_t> && expr < epsilon<real_t>;
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

    using real_t = Real;
    auto vert_to_p0 = tetr_p0 - point;
    auto vert_to_p1 = tetr_p1 - point;
    auto vert_to_p2 = tetr_p2 - point;
    auto vert_to_p3 = tetr_p3 - point;

    std::array<real_t, 5> abs_mixed_prods;
    abs_mixed_prods[0] = std::abs(spt::mixed(vert_to_p0, vert_to_p2, vert_to_p3));
    abs_mixed_prods[1] = std::abs(spt::mixed(vert_to_p0, vert_to_p1, vert_to_p2));
    abs_mixed_prods[2] = std::abs(spt::mixed(vert_to_p0, vert_to_p1, vert_to_p3));
    abs_mixed_prods[3] = std::abs(spt::mixed(vert_to_p1, vert_to_p2, vert_to_p3));
    abs_mixed_prods[4] = std::abs(spt::mixed(tetr_p1 - tetr_p0, tetr_p2 - tetr_p0, tetr_p3 - tetr_p0));

    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4];
}

template <typename Real>
vec<3, Real> closest_segment_point_to_point(
    const vec<3, Real>& point,
    const vec<3, Real>& segm_p0, const vec<3, Real>& segm_p1) {

    using real_t = vec<3>::value_type;
    auto proj = project(point, segm_p0, segm_p1);

    if (helpers::in_rectangle(segm_p0, segm_p1, proj)) {
        return proj;
    } else if (std::array<real_t, 2> sqr_magns
               { (segm_p0 - point).sqr_magnitude(), (segm_p1 - point).sqr_magnitude() };
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

    using real_t = Real;
    std::array<vec<3>, 3> closest_points;
    closest_points[0] = closest_segment_point_to_point(point, trngl_p0, trngl_p1);
    closest_points[1] = closest_segment_point_to_point(point, trngl_p1, trngl_p2);
    closest_points[2] = closest_segment_point_to_point(point, trngl_p2, trngl_p0);

    std::array<real_t, 3> sqrs;
    sqrs[0] = (closest_points[0] - point).sqr_magnitude();
    sqrs[1] = (closest_points[1] - point).sqr_magnitude();
    sqrs[2] = (closest_points[2] - point).sqr_magnitude();

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

    using real_t = Real;
    auto proj = project(point, segm_p0, segm_p1);

    if (helpers::in_rectangle(segm_p0, segm_p1, proj)) {
        return (proj - point).magnitude();
    } else if (std::array<real_t, 2> sqr_magns{ (segm_p0 - point).sqr_magnitude(), (segm_p1 - point).sqr_magnitude() };
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

    using real_t = Real;
    std::array<vec<3>, 3> closest_points;
    closest_points[0] = closest_segment_point_to_point(point, trngl_p0, trngl_p1);
    closest_points[1] = closest_segment_point_to_point(point, trngl_p1, trngl_p2);
    closest_points[2] = closest_segment_point_to_point(point, trngl_p2, trngl_p0);

    std::array<real_t, 3> sqrs;
    sqrs[0] = (closest_points[0] - point).sqr_magnitude();
    sqrs[1] = (closest_points[1] - point).sqr_magnitude();
    sqrs[2] = (closest_points[2] - point).sqr_magnitude();

    return std::min({ sqrs[0], sqrs[1], sqrs[2] });
}

template <typename Real>
vec<3, Real> lines_closest_point(
    const vec<3, Real>& line0_p0, const vec<3, Real>& line0_p1,
    const vec<3, Real>& line1_p0, const vec<3, Real>& line1_p1) {

    using real_t = Real;
    auto u = line0_p1 - line0_p0;
    auto v = line1_p1 - line1_p0;
    auto w = line0_p0 - line1_p0;
    auto a = spt::dot(u, u); // >= 0
    auto b = spt::dot(u, v);
    auto c = spt::dot(v, v); // >= 0
    auto d = spt::dot(u, w);
    auto e = spt::dot(v, w);
    auto determ = helpers::det(a, b, b, c); // >= 0

    real_t sc, tc;
    if (determ < epsilon<real_t>) {
        sc = static_cast<real_t>(0);
        tc = b > c ? d / b : e / c;
    } else {
        auto inv_determ = static_cast<real_t>(1) / determ;
        sc = helpers::det(b, c, d, e) * inv_determ;
        tc = helpers::det(a, b, d, e) * inv_determ;
    }

    return (line0_p0 + u * sc + line1_p0 + v * tc) * static_cast<real_t>(0.5);
}

template <typename Real>
Real lines_distance(
    const vec<3, Real>& line0_p0, const vec<3, Real>& line0_p1,
    const vec<3, Real>& line1_p0, const vec<3, Real>& line1_p1) {

    using real_t = Real;
    auto u = line0_p1 - line0_p0;
    auto v = line1_p1 - line1_p0;
    auto w = line0_p0 - line1_p0;
    auto a = spt::dot(u, u); // >= 0
    auto b = spt::dot(u, v);
    auto c = spt::dot(v, v); // >= 0
    auto d = spt::dot(u, w);
    auto e = spt::dot(v, w);
    auto determ = helpers::det(a, b, b, c); // >= 0

    real_t sc, tc;
    if (determ < epsilon<real_t>) {
        sc = static_cast<real_t>(0);
        tc = b > c ? d / b : e / c;
    } else {
        auto inv_determ = static_cast<real_t>(1) / determ;
        sc = helpers::det(b, c, d, e) * inv_determ;
        tc = helpers::det(a, b, d, e) * inv_determ;
    }

    auto diff_p = w + (u * sc) - (v * tc);
    return diff_p.magnitude();
}

template <typename Real>
Real segments_distance(
    const vec<3, Real>& segm0_p0, const vec<3, Real>& segm0_p1,
    const vec<3, Real>& segm1_p0, const vec<3, Real>& segm1_p1) {

    using real_t = Real;
    auto u = segm0_p1 - segm0_p0;
    auto v = segm1_p1 - segm1_p0;
    auto w = segm0_p0 - segm1_p0;
    auto a = spt::dot(u, u); // >= 0
    auto b = spt::dot(u, v);
    auto c = spt::dot(v, v); // >= 0
    auto d = spt::dot(u, w);
    auto e = spt::dot(v, w);
    auto determ = helpers::det(a, b, b, c); // >= 0
    real_t sc, sn, sd = determ;
    real_t tc, tn, td = determ;

    if (determ < epsilon<real_t>) {
        sn = static_cast<real_t>(0);
        sd = static_cast<real_t>(1);
        tn = e;
        td = c;
    } else {
        sn = helpers::det(b, c, d, e);
        tn = helpers::det(a, b, d, e);

        if (sn < static_cast<real_t>(0)) {
            sn = static_cast<real_t>(0);
            tn = e;
            td = c;
        } else if (sn > sd) {
            sn = sd;
            tn = e + b;
            td = c;
        }
    }

    if (tn < static_cast<real_t>(0)) {
        tn = static_cast<real_t>(0);

        if (-d < static_cast<real_t>(0)) {
            sn = static_cast<real_t>(0);
        } else if (-d > a) {
            sn = sd;
        } else {
            sn = -d;
            sd = a;
        }
    } else if (tn > td) {
        tn = td;

        if (b - d < static_cast<real_t>(0)) {
            sn = static_cast<real_t>(0);
        } else if (b - d > a) {
            sn = sd;
        } else {
            sn = b - d;
            sd = a;
        }
    }

    sc = std::abs(sn) < epsilon<real_t> ? static_cast<real_t>(0) : sn / sd;
    tc = std::abs(tn) < epsilon<real_t> ? static_cast<real_t>(0) : tn / td;

    auto diff_p = w + (u * sc) - (v * tc);

    return diff_p.magnitude();
}

// CPA - Closest Point of Approach.
template <typename Real>
Real cpa_time(
    const vec<3, Real>& start0, const vec<3, Real>& vel0,
    const vec<3, Real>& start1, const vec<3, Real>& vel1) {

    using real_t = vec<3>::value_type;
    auto dv = vel0 - vel1;
    auto dv2 = spt::dot(dv, dv);
    if (dv2 < epsilon<real_t>)
        return static_cast<real_t>(0);

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

} // namespace spt::algs
