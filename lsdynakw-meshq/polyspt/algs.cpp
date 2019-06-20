// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

// With help of GeomAlgorithms.com website.

#include "algs.h"

// TODO: remove defines
#define EPS static_cast<vec3::real_type>(1e-6)

#define DET(a, b, c, d) \
        ((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define IN_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0.x[0], corner1.x[0], point.x[0]) && \
         BETWEEN(corner0.x[1], corner1.x[1], point.x[1]) && \
         BETWEEN(corner0.x[2], corner1.x[2], point.x[2]))


using spt::algs::vec3;


vec3 spt::algs::project(
        const vec3& point,
        const vec3& line_p0, const vec3& line_p1)
{
    return line_p0 + (point - line_p0).project(line_p1 - line_p0);
}


bool spt::algs::project(
    vec3& out,
    const vec3& point,
    const vec3& segm_p0, const vec3& segm_p1)
{
    vec3 res = segm_p0 + (point - segm_p0).project(segm_p1 - segm_p0);

    if (IN_RECTANGLE(
            segm_p0,
            segm_p1,
            res))
    {
        out = res;
        return true;
    }

    return false;
}


vec3 spt::algs::project(
        const vec3& point,
        const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2)
{
    return plane_p0 + (point - plane_p0).project(plane_p1 - plane_p0, plane_p2 - plane_p0);
}




bool spt::algs::does_ray_intersect_plane(
        const vec3& dir,
        const vec3& pl_p0, const vec3& pl_p1, const vec3& pl_p2)
{
    std::array<vec3, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    vec3::real_type det = vec3::dot(edges[0], pvec);
    return det > static_cast<vec3::real_type>(1e-6) || det < static_cast<vec3::real_type>(-1e-6);
}


bool spt::algs::ray_intersect_plane(
    vec3& out_intersectPoint,
    const vec3& origin, const vec3& dir,
    const vec3& pl_p0, const vec3& pl_p1, const vec3& pl_p2)
{
    std::array<vec3, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    vec3::real_type det = vec3::dot(edges[0], pvec);
    if (det < static_cast<vec3::real_type>(1e-6) && det > static_cast<vec3::real_type>(-1e-6))
        return false;

    vec3 tvec = origin - pl_p0;
    vec3 qvec = vec3::cross(tvec, edges[0]);

    vec3::real_type t = vec3::dot(edges[1], qvec) / det;
    out_intersectPoint = origin + dir * t;
    return t > static_cast<vec3::real_type>(0.0);
}


bool spt::algs::does_ray_intersect_triangle(
    const vec3& origin, const vec3& dir,
    const vec3& tr_p0, const vec3& tr_p1, const vec3& tr_p2)
{
    std::array<vec3, 2> edges = { tr_p1 - tr_p0, tr_p2 - tr_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    vec3::real_type det = vec3::dot(edges[0], pvec);
    if (det < static_cast<vec3::real_type>(1e-6) && det > static_cast<vec3::real_type>(-1e-6))
        return false;

    vec3::real_type inv_det = static_cast<vec3::real_type>(1.0) / det;

    vec3 tvec = origin - tr_p0;
    vec3::real_type u = vec3::dot(tvec, pvec) * inv_det;
    if (u < static_cast<vec3::real_type>(0.0) || u > static_cast<vec3::real_type>(1.0))
        return false;

    vec3 qvec = vec3::cross(tvec, edges[0]);
    vec3::real_type v = vec3::dot(dir, qvec) * inv_det;
    if (v < static_cast<vec3::real_type>(0.0) || u + v > static_cast<vec3::real_type>(1.0))
        return false;

    vec3::real_type t = vec3::dot(edges[1], qvec) * inv_det;
    return t > static_cast<vec3::real_type>(0.0);
}


vec3 spt::algs::lines_closest_point(const vec3& line0_p0, const vec3& line0_p1, const vec3& line1_p0, const vec3& line1_p1)
{
    vec3 u = line0_p1 - line0_p0;
    vec3 v = line1_p1 - line1_p0;
    vec3 w = line0_p0 - line1_p0;
    vec3::real_type a = vec3::dot(u, u); // >= 0
    vec3::real_type b = vec3::dot(u, v);
    vec3::real_type c = vec3::dot(v, v); // >= 0
    vec3::real_type d = vec3::dot(u, w);
    vec3::real_type e = vec3::dot(v, w);
    vec3::real_type determ = DET(a, b, b, c); // >= 0
    vec3::real_type sc, tc;

    if (determ < static_cast<vec3::real_type>(1e-6))
    {
        sc = static_cast<vec3::real_type>(0.0);
        tc = b > c ? d / b : e / c;
    }
    else
    {
        vec3::real_type inv_determ = static_cast<vec3::real_type>(1.0) / determ;
        sc = DET(b, c, d, e) * inv_determ;
        tc = DET(a, b, d, e) * inv_determ;
    }

    return (line0_p0 + u * sc + line1_p0 + v * tc) * static_cast<vec3::real_type>(0.5);
}


bool spt::algs::line_intersect_plane(
    vec3& out_intersectPoint,
    const vec3& line_point, const vec3& line_dir,
    const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2)
{
    std::array<vec3, 2> edges = { plane_p1 - plane_p0, plane_p2 - plane_p0 };

    vec3 pvec = vec3::cross(line_dir, edges[1]);
    vec3::real_type det = vec3::dot(edges[0], pvec);
    if (det < static_cast<vec3::real_type>(1e-6) && det > static_cast<vec3::real_type>(-1e-6))
        return false;

    vec3 tvec = line_point - plane_p0;
    vec3 qvec = vec3::cross(tvec, edges[0]);

    vec3::real_type t = vec3::dot(edges[1], qvec) / det;
    out_intersectPoint = line_point + line_dir * t;
    return true;
}


vec3 spt::algs::line_intersect_plane(
    const vec3& line_point, const vec3& line_dir,
    const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2)
{
    std::array<vec3, 2> edges = { plane_p1 - plane_p0, plane_p2 - plane_p0 };

    vec3 pvec = vec3::cross(line_dir, edges[1]);
    vec3 tvec = line_point - plane_p0;
    vec3 qvec = vec3::cross(tvec, edges[0]);

    vec3::real_type t = vec3::dot(edges[1], qvec) / vec3::dot(edges[0], pvec);
    return line_point + line_dir * t;
}


bool spt::algs::does_segment_intersect_triangle(
    const vec3& segm_p0, const vec3& segm_p1,
    const vec3& tr_p0, const vec3& tr_p1, const vec3& tr_p2)
{
    const vec3 dir = segm_p1 - segm_p0;

    std::array<vec3, 2> edges = { tr_p1 - tr_p0, tr_p2 - tr_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    vec3::real_type det = vec3::dot(edges[0], pvec);
    if (det < static_cast<vec3::real_type>(1e-6) && det > static_cast<vec3::real_type>(-1e-6))
        return false;

    vec3::real_type inv_det = static_cast<vec3::real_type>(1.0) / det;

    vec3 tvec = segm_p0 - tr_p0;
    vec3::real_type u = vec3::dot(tvec, pvec) * inv_det;
    if (u < static_cast<vec3::real_type>(0.0) || u > static_cast<vec3::real_type>(1.0))
        return false;

    vec3 qvec = vec3::cross(tvec, edges[0]);
    vec3::real_type v = vec3::dot(dir, qvec) * inv_det;
    if (v < static_cast<vec3::real_type>(0.0) || u + v > static_cast<vec3::real_type>(1.0))
        return false;

    vec3::real_type t = vec3::dot(edges[1], qvec) * inv_det;
    return t < static_cast<vec3::real_type>(1.0) && t > static_cast<vec3::real_type>(0.0);
}


bool spt::algs::segment_intersect_plane(
    vec3& out_intersectPoint,
    const vec3& p0, const vec3& p1,
    const vec3& pl_p0, const vec3& pl_p1, const vec3& pl_p2)
{
    const vec3 dir = p1 - p0;

    std::array<vec3, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    vec3::real_type det = vec3::dot(edges[0], pvec);
    if (det < static_cast<vec3::real_type>(1e-6) && det > static_cast<vec3::real_type>(-1e-6))
        return false;

    vec3 tvec = p0 - pl_p0;
    vec3 qvec = vec3::cross(tvec, edges[0]);

    vec3::real_type t = vec3::dot(edges[1], qvec) / det;
    out_intersectPoint = p0 + dir * t;
    return t < static_cast<vec3::real_type>(1.0) && t > static_cast<vec3::real_type>(0.0);
}




bool spt::algs::is_point_on_triangle(
        const vec3& point,
        const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    vec3::real_type s0 = vec3::cross(trngl_p0 - point,    trngl_p1 - point   ).magnitude();
    vec3::real_type s1 = vec3::cross(trngl_p0 - point,    trngl_p2 - point   ).magnitude();
    vec3::real_type s2 = vec3::cross(trngl_p1 - point,    trngl_p2 - point   ).magnitude();
    vec3::real_type s  = vec3::cross(trngl_p0 - trngl_p2, trngl_p1 - trngl_p2).magnitude();

    vec3::real_type expr = s - s0 - s1 - s2;
    return expr > static_cast<vec3::real_type>(-1e-6) && expr < static_cast<vec3::real_type>(1e-6);
}


bool spt::algs::is_point_on_triangle(
        const vec3& point,
        const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2,
        vec3::real_type max_sqrs_sum)
{
    if (sqrs_sum(point, trngl_p0, trngl_p1, trngl_p2) > max_sqrs_sum)
        return false;

    return is_point_on_triangle(point, trngl_p0, trngl_p1, trngl_p2);
}




vec3::real_type spt::algs::lines_distance(
        const vec3& line0_p0, const vec3& line0_p1,
        const vec3& line1_p0, const vec3& line1_p1)
{
    vec3 u = line0_p1 - line0_p0;
    vec3 v = line1_p1 - line1_p0;
    vec3 w = line0_p0 - line1_p0;
    vec3::real_type a = vec3::dot(u, u); // >= 0
    vec3::real_type b = vec3::dot(u, v);
    vec3::real_type c = vec3::dot(v, v); // >= 0
    vec3::real_type d = vec3::dot(u, w);
    vec3::real_type e = vec3::dot(v, w);
    vec3::real_type determ = DET(a, b, b, c); // >= 0
    vec3::real_type sc, tc;

    if (determ < static_cast<vec3::real_type>(1e-6))
    {
        sc = 0.0;
        tc = b > c ? d / b : e / c;
    }
    else
    {
        vec3::real_type inv_determ = static_cast<vec3::real_type>(1.0) / determ;
        sc = DET(b, c, d, e) * inv_determ;
        tc = DET(a, b, d, e) * inv_determ;
    }

    vec3 diff_p = w + (u * sc) - (v * tc);

    return diff_p.magnitude();
}


vec3::real_type spt::algs::segments_distance(
    const vec3& segm0_p0, const vec3& segm0_p1,
    const vec3& segm1_p0, const vec3& segm1_p1)
{
    vec3 u = segm0_p1 - segm0_p0;
    vec3 v = segm1_p1 - segm1_p0;
    vec3 w = segm0_p0 - segm1_p0;
    vec3::real_type a = vec3::dot(u, u); // >= 0
    vec3::real_type b = vec3::dot(u, v);
    vec3::real_type c = vec3::dot(v, v); // >= 0
    vec3::real_type d = vec3::dot(u, w);
    vec3::real_type e = vec3::dot(v, w);
    vec3::real_type determ = DET(a, b, b, c); // >= 0
    vec3::real_type sc, sn, sd = determ;
    vec3::real_type tc, tn, td = determ;

    if (determ < static_cast<vec3::real_type>(1e-6))
    {
        sn = static_cast<vec3::real_type>(0.0);
        sd = static_cast<vec3::real_type>(1.0);
        tn = e;
        td = c;
    }
    else
    {
        sn = DET(b, c, d, e);
        tn = DET(a, b, d, e);

        if (sn < static_cast<vec3::real_type>(0.0))
        {
            sn = 0.0;
            tn = e;
            td = c;
        }
        else if (sn > sd)
        {
            sn = sd;
            tn = e + b;
            td = c;
        }
    }

    if (tn < static_cast<vec3::real_type>(0.0))
    {
        tn = static_cast<vec3::real_type>(0.0);

        if (-d < static_cast<vec3::real_type>(0.0))
        {
            sn = static_cast<vec3::real_type>(0.0);
        }
        else if (-d > a)
        {
            sn = sd;
        }
        else
        {
            sn = -d;
            sd = a;
        }
    }
    else if (tn > td)
    {
        tn = td;

        if (b - d < static_cast<vec3::real_type>(0.0))
        {
            sn = static_cast<vec3::real_type>(0.0);
        }
        else if (b - d > a)
        {
            sn = sd;
        }
        else
        {
            sn = b - d;
            sd = a;
        }
    }

    sc = std::abs(sn) < static_cast<vec3::real_type>(1e-6) ? static_cast<vec3::real_type>(0.0) : sn / sd;
    tc = std::abs(tn) < static_cast<vec3::real_type>(1e-6) ? static_cast<vec3::real_type>(0.0) : tn / td;

    vec3 diff_p = w + (u * sc) - (v * tc);

    return diff_p.magnitude();
}




vec3::real_type spt::algs::cpa_time(
    const vec3& start0, const vec3& vel0,
    const vec3& start1, const vec3& vel1)
{
    vec3 dv = vel0 - vel1;
    vec3::real_type dv2 = vec3::dot(dv, dv);
    if (dv2 < static_cast<vec3::real_type>(1e-6))
        return static_cast<vec3::real_type>(0.0);

    vec3 w0 = start0 - start1;
    return -vec3::dot(w0, dv) / dv2;
}



vec3::real_type spt::algs::cpa_distance(
    const vec3& start0, const vec3& vel0,
    const vec3& start1, const vec3& vel1)
{
    vec3::real_type time = cpa_time(start0, vel0, start1, vel1);
    vec3 p0 = start0 + vel0 * time;
    vec3 p1 = start1 + vel1 * time;
    return (p1 - p0).magnitude();
}




vec3::real_type spt::algs::distance_point_to_line(const vec3& point, const vec3& line_p0, const vec3& line_p1)
{
    return (project(point, line_p0, line_p1) - point).magnitude();
}


vec3::real_type spt::algs::distance_point_to_segment(const vec3& point, const vec3& segm_p0, const vec3& segm_p1)
{
    vec3 proj = project(point, segm_p0, segm_p1);

    if (IN_RECTANGLE(segm_p0, segm_p1, proj))
    {
        return (proj - point).magnitude();
    }
    else if (vec3::real_type sqr_magns[2] { (segm_p0 - point).sqrMagnitude(), (segm_p1 - point).sqrMagnitude() };
             sqr_magns[0] < sqr_magns[1])
    {
        return std::sqrt(sqr_magns[0]);
    }
    else
    {
        return std::sqrt(sqr_magns[1]);
    }
}


vec3::real_type spt::algs::distance_point_to_triangle_on_plane(
        const vec3& point, const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    std::array<vec3, 3> closest_points;
    closest_points[0] = closest_segment_point_to_point(point, trngl_p0, trngl_p1);
    closest_points[1] = closest_segment_point_to_point(point, trngl_p1, trngl_p2);
    closest_points[2] = closest_segment_point_to_point(point, trngl_p2, trngl_p0);

    std::array<vec3::real_type, 3> sqrs;
    sqrs[0] = (closest_points[0] - point).sqrMagnitude();
    sqrs[1] = (closest_points[1] - point).sqrMagnitude();
    sqrs[2] = (closest_points[2] - point).sqrMagnitude();

    return std::min({ sqrs[0], sqrs[1], sqrs[2] });
}


bool spt::algs::does_triangle_intersect_sphere(const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2,
                                           const vec3& center,   vec3::real_type radius)
{
    vec3 proj = project(center, trngl_p0, trngl_p1, trngl_p2);
    if ((proj - center).sqrMagnitude() > radius * radius)
        return false;

    if (is_point_on_triangle(proj, trngl_p0, trngl_p1, trngl_p2, max_sqrs_sum(trngl_p0, trngl_p1, trngl_p2)))
        return true;

    vec3 closest = closest_triangle_point_to_point_on_plane(proj, trngl_p0, trngl_p1, trngl_p2);
    return (closest - center).sqrMagnitude() <= radius * radius;
}


vec3 spt::algs::closest_segment_point_to_point(const vec3& point, const vec3& segm_p0, const vec3& segm_p1)
{
    vec3 proj = project(point, segm_p0, segm_p1);

    if (IN_RECTANGLE(segm_p0, segm_p1, proj))
    {
        return proj;
    }
    else if (vec3::real_type sqr_magns[2] { (segm_p0 - point).sqrMagnitude(), (segm_p1 - point).sqrMagnitude() };
             sqr_magns[0] < sqr_magns[1])
    {
        return segm_p0;
    }
    else
    {
        return segm_p1;
    }
}


vec3::real_type spt::algs::max_sqrs_sum(const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    std::array<vec3::real_type, 3> sqrs;
    sqrs[0] = (trngl_p1 - trngl_p0).sqrMagnitude();
    sqrs[1] = (trngl_p2 - trngl_p1).sqrMagnitude();
    sqrs[2] = (trngl_p0 - trngl_p2).sqrMagnitude();

    std::size_t max_inds[2];
    if (sqrs[0] < sqrs[1])
    {
        max_inds[0] = 1;
        if (sqrs[0] < sqrs[2])
            max_inds[1] = 2;
        else
            max_inds[1] = 0;
    }
    else
    {
        max_inds[0] = 0;
        if (sqrs[1] < sqrs[2])
            max_inds[1] = 2;
        else
            max_inds[1] = 1;
    }

    return sqrs[max_inds[0]] + sqrs[max_inds[1]];
}


vec3::real_type spt::algs::sqrs_sum(const vec3& point, const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    std::array<vec3::real_type, 3> sqrs;
    sqrs[0] = (trngl_p0 - point).sqrMagnitude();
    sqrs[1] = (trngl_p1 - point).sqrMagnitude();
    sqrs[2] = (trngl_p2 - point).sqrMagnitude();
    return sqrs[0] + sqrs[1] + sqrs[2];
}


vec3 spt::algs::closest_triangle_point_to_point_on_plane(
        const vec3& point, const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    std::array<vec3, 3> closest_points;
    closest_points[0] = closest_segment_point_to_point(point, trngl_p0, trngl_p1);
    closest_points[1] = closest_segment_point_to_point(point, trngl_p1, trngl_p2);
    closest_points[2] = closest_segment_point_to_point(point, trngl_p2, trngl_p0);

    std::array<vec3::real_type, 3> sqrs;
    sqrs[0] = (closest_points[0] - point).sqrMagnitude();
    sqrs[1] = (closest_points[1] - point).sqrMagnitude();
    sqrs[2] = (closest_points[2] - point).sqrMagnitude();

    std::size_t min_i;
    if (sqrs[0] < sqrs[1])
    {
        if (sqrs[0] < sqrs[2])
            min_i = 0;
        else
            min_i = 2;
    }
    else
    {
        if (sqrs[1] < sqrs[2])
            min_i = 1;
        else
            min_i = 2;
    }

    return closest_points[min_i];
}


bool spt::algs::is_point_in_tetrahedron(
        const vec3& point, const vec3& tetr_p0, const vec3& tetr_p1, const vec3& tetr_p2, const vec3& tetr_p3)
{
    vec3 vert_to_p0 = tetr_p0 - point;
    vec3 vert_to_p1 = tetr_p1 - point;
    vec3 vert_to_p2 = tetr_p2 - point;
    vec3 vert_to_p3 = tetr_p3 - point;

    std::array<vec3::real_type, 5> abs_mixed_prods;
    abs_mixed_prods[0] = std::abs(vec3::mixed(vert_to_p0, vert_to_p2, vert_to_p3));
    abs_mixed_prods[1] = std::abs(vec3::mixed(vert_to_p0, vert_to_p1, vert_to_p2));
    abs_mixed_prods[2] = std::abs(vec3::mixed(vert_to_p0, vert_to_p1, vert_to_p3));
    abs_mixed_prods[3] = std::abs(vec3::mixed(vert_to_p1, vert_to_p2, vert_to_p3));
    abs_mixed_prods[4] = std::abs(vec3::mixed(tetr_p1 - tetr_p0, tetr_p2 - tetr_p0, tetr_p3 - tetr_p0));

    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4];
}
