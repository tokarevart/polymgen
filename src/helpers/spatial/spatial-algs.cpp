// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

// With help of GeomAlgorithms.com website.

#include "algs.h"


#define EPS static_cast<real_t>(1e-6)

#define DET(a, b, c, d) \
        ((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define IN_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0.x[0], corner1.x[0], point.x[0]) && \
         BETWEEN(corner0.x[1], corner1.x[1], point.x[1]) && \
         BETWEEN(corner0.x[2], corner1.x[2], point.x[2]))


using spt::vec3;


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




bool spt::algs::doesRayIntersectPlane(
        const vec3& dir,
        const vec3& pl_p0, const vec3& pl_p1, const vec3& pl_p2)
{
    std::array<vec3, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    real_t det = vec3::dot(edges[0], pvec);
    return det > static_cast<real_t>(1e-6) || det < static_cast<real_t>(-1e-6);
}


bool spt::algs::rayIntersectPlane(
    vec3& out_intersectPoint,
    const vec3& origin, const vec3& dir,
    const vec3& pl_p0, const vec3& pl_p1, const vec3& pl_p2)
{
    std::array<vec3, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    real_t det = vec3::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    vec3 tvec = origin - pl_p0;
    vec3 qvec = vec3::cross(tvec, edges[0]);

    real_t t = vec3::dot(edges[1], qvec) / det;
    out_intersectPoint = origin + dir * t;
    return t > static_cast<real_t>(0.0);
}


bool spt::algs::doesRayIntersectTriangle(
    const vec3& origin, const vec3& dir,
    const vec3& tr_p0, const vec3& tr_p1, const vec3& tr_p2)
{
    std::array<vec3, 2> edges = { tr_p1 - tr_p0, tr_p2 - tr_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    real_t det = vec3::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    real_t inv_det = static_cast<real_t>(1.0) / det;

    vec3 tvec = origin - tr_p0;
    real_t u = vec3::dot(tvec, pvec) * inv_det;
    if (u < static_cast<real_t>(0.0) || u > static_cast<real_t>(1.0))
        return false;

    vec3 qvec = vec3::cross(tvec, edges[0]);
    real_t v = vec3::dot(dir, qvec) * inv_det;
    if (v < static_cast<real_t>(0.0) || u + v > static_cast<real_t>(1.0))
        return false;

    real_t t = vec3::dot(edges[1], qvec) * inv_det;
    return t > static_cast<real_t>(0.0);
}


vec3 spt::algs::linesClosestPoint(const vec3& line0_p0, const vec3& line0_p1, const vec3& line1_p0, const vec3& line1_p1)
{
    vec3 u = line0_p1 - line0_p0;
    vec3 v = line1_p1 - line1_p0;
    vec3 w = line0_p0 - line1_p0;
    real_t a = vec3::dot(u, u); // >= 0
    real_t b = vec3::dot(u, v);
    real_t c = vec3::dot(v, v); // >= 0
    real_t d = vec3::dot(u, w);
    real_t e = vec3::dot(v, w);
    real_t determ = DET(a, b, b, c); // >= 0
    real_t sc, tc;

    if (determ < static_cast<real_t>(1e-6))
    {
        sc = static_cast<real_t>(0.0);
        tc = b > c ? d / b : e / c;
    }
    else
    {
        real_t inv_determ = static_cast<real_t>(1.0) / determ;
        sc = DET(b, c, d, e) * inv_determ;
        tc = DET(a, b, d, e) * inv_determ;
    }

    return (line0_p0 + u * sc + line1_p0 + v * tc) * static_cast<real_t>(0.5);
}


bool spt::algs::lineIntersectPlane(
    vec3& out_intersectPoint,
    const vec3& line_point, const vec3& line_dir,
    const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2)
{
    std::array<vec3, 2> edges = { plane_p1 - plane_p0, plane_p2 - plane_p0 };

    vec3 pvec = vec3::cross(line_dir, edges[1]);
    real_t det = vec3::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    vec3 tvec = line_point - plane_p0;
    vec3 qvec = vec3::cross(tvec, edges[0]);

    real_t t = vec3::dot(edges[1], qvec) / det;
    out_intersectPoint = line_point + line_dir * t;
    return true;
}


vec3 spt::algs::lineIntersectPlane(
    const vec3& line_point, const vec3& line_dir,
    const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2)
{
    std::array<vec3, 2> edges = { plane_p1 - plane_p0, plane_p2 - plane_p0 };

    vec3 pvec = vec3::cross(line_dir, edges[1]);
    vec3 tvec = line_point - plane_p0;
    vec3 qvec = vec3::cross(tvec, edges[0]);

    real_t t = vec3::dot(edges[1], qvec) / vec3::dot(edges[0], pvec);
    return line_point + line_dir * t;
}


bool spt::algs::doesSegmentIntersectTriangle(
    const vec3& segm_p0, const vec3& segm_p1,
    const vec3& tr_p0, const vec3& tr_p1, const vec3& tr_p2)
{
    const vec3 dir = segm_p1 - segm_p0;

    std::array<vec3, 2> edges = { tr_p1 - tr_p0, tr_p2 - tr_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    real_t det = vec3::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    real_t inv_det = static_cast<real_t>(1.0) / det;

    vec3 tvec = segm_p0 - tr_p0;
    real_t u = vec3::dot(tvec, pvec) * inv_det;
    if (u < static_cast<real_t>(0.0) || u > static_cast<real_t>(1.0))
        return false;

    vec3 qvec = vec3::cross(tvec, edges[0]);
    real_t v = vec3::dot(dir, qvec) * inv_det;
    if (v < static_cast<real_t>(0.0) || u + v > static_cast<real_t>(1.0))
        return false;

    real_t t = vec3::dot(edges[1], qvec) * inv_det;
    return t < static_cast<real_t>(1.0) && t > static_cast<real_t>(0.0);
}


bool spt::algs::segmentIntersectPlane(
    vec3& out_intersectPoint,
    const vec3& p0, const vec3& p1,
    const vec3& pl_p0, const vec3& pl_p1, const vec3& pl_p2)
{
    const vec3 dir = p1 - p0;

    std::array<vec3, 2> edges = { pl_p1 - pl_p0, pl_p2 - pl_p0 };

    vec3 pvec = vec3::cross(dir, edges[1]);
    real_t det = vec3::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    vec3 tvec = p0 - pl_p0;
    vec3 qvec = vec3::cross(tvec, edges[0]);

    real_t t = vec3::dot(edges[1], qvec) / det;
    out_intersectPoint = p0 + dir * t;
    return t < static_cast<real_t>(1.0) && t > static_cast<real_t>(0.0);
}




bool spt::algs::isPointOnTriangle(
        const vec3& point,
        const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    real_t s0 = vec3::cross(trngl_p0 - point,    trngl_p1 - point   ).magnitude();
    real_t s1 = vec3::cross(trngl_p0 - point,    trngl_p2 - point   ).magnitude();
    real_t s2 = vec3::cross(trngl_p1 - point,    trngl_p2 - point   ).magnitude();
    real_t s  = vec3::cross(trngl_p0 - trngl_p2, trngl_p1 - trngl_p2).magnitude();

    real_t expr = s - s0 - s1 - s2;
    return expr > static_cast<real_t>(-1e-6) && expr < static_cast<real_t>(1e-6);
}


bool spt::algs::isPointOnTriangle(
        const vec3& point,
        const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2,
        real_t max_sqrs_sum)
{
    if (computeSqrsSum(point, trngl_p0, trngl_p1, trngl_p2) > max_sqrs_sum)
        return false;

    return isPointOnTriangle(point, trngl_p0, trngl_p1, trngl_p2);
}




real_t spt::algs::linesDistance(
        const vec3& line0_p0, const vec3& line0_p1,
        const vec3& line1_p0, const vec3& line1_p1)
{
    vec3 u = line0_p1 - line0_p0;
    vec3 v = line1_p1 - line1_p0;
    vec3 w = line0_p0 - line1_p0;
    real_t a = vec3::dot(u, u); // >= 0
    real_t b = vec3::dot(u, v);
    real_t c = vec3::dot(v, v); // >= 0
    real_t d = vec3::dot(u, w);
    real_t e = vec3::dot(v, w);
    real_t determ = DET(a, b, b, c); // >= 0
    real_t sc, tc;

    if (determ < static_cast<real_t>(1e-6))
    {
        sc = 0.0;
        tc = b > c ? d / b : e / c;
    }
    else
    {
        real_t inv_determ = static_cast<real_t>(1.0) / determ;
        sc = DET(b, c, d, e) * inv_determ;
        tc = DET(a, b, d, e) * inv_determ;
    }

    vec3 diff_p = w + (u * sc) - (v * tc);

    return diff_p.magnitude();
}


real_t spt::algs::segmentsDistance(
    const vec3& segm0_p0, const vec3& segm0_p1,
    const vec3& segm1_p0, const vec3& segm1_p1)
{
    vec3 u = segm0_p1 - segm0_p0;
    vec3 v = segm1_p1 - segm1_p0;
    vec3 w = segm0_p0 - segm1_p0;
    real_t a = vec3::dot(u, u); // >= 0
    real_t b = vec3::dot(u, v);
    real_t c = vec3::dot(v, v); // >= 0
    real_t d = vec3::dot(u, w);
    real_t e = vec3::dot(v, w);
    real_t determ = DET(a, b, b, c); // >= 0
    real_t sc, sn, sd = determ;
    real_t tc, tn, td = determ;

    if (determ < static_cast<real_t>(1e-6))
    {
        sn = static_cast<real_t>(0.0);
        sd = static_cast<real_t>(1.0);
        tn = e;
        td = c;
    }
    else
    {
        sn = DET(b, c, d, e);
        tn = DET(a, b, d, e);

        if (sn < static_cast<real_t>(0.0))
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

    if (tn < static_cast<real_t>(0.0))
    {
        tn = static_cast<real_t>(0.0);

        if (-d < static_cast<real_t>(0.0))
        {
            sn = static_cast<real_t>(0.0);
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

        if (b - d < static_cast<real_t>(0.0))
        {
            sn = static_cast<real_t>(0.0);
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

    sc = std::abs(sn) < static_cast<real_t>(1e-6) ? static_cast<real_t>(0.0) : sn / sd;
    tc = std::abs(tn) < static_cast<real_t>(1e-6) ? static_cast<real_t>(0.0) : tn / td;

    vec3 diff_p = w + (u * sc) - (v * tc);

    return diff_p.magnitude();
}




real_t spt::algs::cpaTime(
    const vec3& start0, const vec3& vel0,
    const vec3& start1, const vec3& vel1)
{
    vec3 dv = vel0 - vel1;
    real_t dv2 = vec3::dot(dv, dv);
    if (dv2 < static_cast<real_t>(1e-6))
        return static_cast<real_t>(0.0);

    vec3 w0 = start0 - start1;
    return -vec3::dot(w0, dv) / dv2;
}



real_t spt::algs::cpaDistance(
    const vec3& start0, const vec3& vel0,
    const vec3& start1, const vec3& vel1)
{
    real_t cpa_time = cpaTime(start0, vel0, start1, vel1);
    vec3 p0 = start0 + vel0 * cpa_time;
    vec3 p1 = start1 + vel1 * cpa_time;
    return (p1 - p0).magnitude();
}




real_t spt::algs::distancePointToLine(const vec3& point, const vec3& line_p0, const vec3& line_p1)
{
    return (project(point, line_p0, line_p1) - point).magnitude();
}


real_t spt::algs::distancePointToSegment(const vec3& point, const vec3& segm_p0, const vec3& segm_p1)
{
    vec3 proj = project(point, segm_p0, segm_p1);

    if (IN_RECTANGLE(segm_p0, segm_p1, proj))
    {
        return (proj - point).magnitude();
    }
    else if (real_t sqr_magns[2] { (segm_p0 - point).sqrMagnitude(), (segm_p1 - point).sqrMagnitude() };
             sqr_magns[0] < sqr_magns[1])
    {
        return std::sqrt(sqr_magns[0]);
    }
    else
    {
        return std::sqrt(sqr_magns[1]);
    }
}


real_t spt::algs::distancePointToTriangleOnPlane(
        const vec3& point, const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    std::array<vec3, 3> closest_points;
    closest_points[0] = closestSegmentPointToPoint(point, trngl_p0, trngl_p1);
    closest_points[1] = closestSegmentPointToPoint(point, trngl_p1, trngl_p2);
    closest_points[2] = closestSegmentPointToPoint(point, trngl_p2, trngl_p0);

    std::array<real_t, 3> sqrs;
    sqrs[0] = (closest_points[0] - point).sqrMagnitude();
    sqrs[1] = (closest_points[1] - point).sqrMagnitude();
    sqrs[2] = (closest_points[2] - point).sqrMagnitude();

    return std::min({ sqrs[0], sqrs[1], sqrs[2] });
}


bool spt::algs::doesTriangleIntersectSphere(const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2,
                                           const vec3& center,   real_t radius)
{
    vec3 proj = project(center, trngl_p0, trngl_p1, trngl_p2);
    if ((proj - center).sqrMagnitude() > radius * radius)
        return false;

    if (isPointOnTriangle(proj, trngl_p0, trngl_p1, trngl_p2, computeMaxSqrsSum(trngl_p0, trngl_p1, trngl_p2)))
        return true;

    vec3 closest = closestTrianglePointToPointOnPlane(proj, trngl_p0, trngl_p1, trngl_p2);
    return (closest - center).sqrMagnitude() <= radius * radius;
}


vec3 spt::algs::closestSegmentPointToPoint(const vec3& point, const vec3& segm_p0, const vec3& segm_p1)
{
    vec3 proj = project(point, segm_p0, segm_p1);

    if (IN_RECTANGLE(segm_p0, segm_p1, proj))
    {
        return proj;
    }
    else if (real_t sqr_magns[2] { (segm_p0 - point).sqrMagnitude(), (segm_p1 - point).sqrMagnitude() };
             sqr_magns[0] < sqr_magns[1])
    {
        return segm_p0;
    }
    else
    {
        return segm_p1;
    }
}


real_t spt::algs::computeMaxSqrsSum(const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    std::array<real_t, 3> sqrs;
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


real_t spt::algs::computeSqrsSum(const vec3& point, const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    std::array<real_t, 3> sqrs;
    sqrs[0] = (trngl_p0 - point).sqrMagnitude();
    sqrs[1] = (trngl_p1 - point).sqrMagnitude();
    sqrs[2] = (trngl_p2 - point).sqrMagnitude();
    return sqrs[0] + sqrs[1] + sqrs[2];
}


vec3 spt::algs::closestTrianglePointToPointOnPlane(
        const vec3& point, const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2)
{
    std::array<vec3, 3> closest_points;
    closest_points[0] = closestSegmentPointToPoint(point, trngl_p0, trngl_p1);
    closest_points[1] = closestSegmentPointToPoint(point, trngl_p1, trngl_p2);
    closest_points[2] = closestSegmentPointToPoint(point, trngl_p2, trngl_p0);

    std::array<real_t, 3> sqrs;
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


bool spt::algs::isPointInTetrahedron(
        const vec3& point, const vec3& tetr_p0, const vec3& tetr_p1, const vec3& tetr_p2, const vec3& tetr_p3)
{
    vec3 vert_to_p0 = tetr_p0 - point;
    vec3 vert_to_p1 = tetr_p1 - point;
    vec3 vert_to_p2 = tetr_p2 - point;
    vec3 vert_to_p3 = tetr_p3 - point;

    std::array<real_t, 5> abs_mixed_prods;
    abs_mixed_prods[0] = std::abs(vec3::mixed(vert_to_p0, vert_to_p2, vert_to_p3));
    abs_mixed_prods[1] = std::abs(vec3::mixed(vert_to_p0, vert_to_p1, vert_to_p2));
    abs_mixed_prods[2] = std::abs(vec3::mixed(vert_to_p0, vert_to_p1, vert_to_p3));
    abs_mixed_prods[3] = std::abs(vec3::mixed(vert_to_p1, vert_to_p2, vert_to_p3));
    abs_mixed_prods[4] = std::abs(vec3::mixed(tetr_p1 - tetr_p0, tetr_p2 - tetr_p0, tetr_p3 - tetr_p0));

    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4];
}
