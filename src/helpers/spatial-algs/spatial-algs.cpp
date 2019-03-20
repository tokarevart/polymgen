// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

// With help of GeomAlgorithms.com website.

#include "helpers/spatial-algs/spatial-algs.h"


#define EPS static_cast<real_t>(1e-6)

#define DET(a, b, c, d) \
        ((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0.coors[0], corner1.coors[0], point.coors[0]) && \
         BETWEEN(corner0.coors[1], corner1.coors[1], point.coors[1]) && \
         BETWEEN(corner0.coors[2], corner1.coors[2], point.coors[2]))





Vec spatalgs::project(
        const Vec& point,
        const Vec& line_p0, const Vec& line_p1)
{
    return line_p0 + (point - line_p0).project(line_p1 - line_p0);
}


bool spatalgs::project(
    Vec& out,
    const Vec& point,
    const Vec& segm_p0, const Vec& segm_p1)
{
    Vec res = segm_p0 + (point - segm_p0).project(segm_p1 - segm_p0);

    if (INSIDE_RECTANGLE(
            segm_p0,
            segm_p1,
            res))
    {
        out = res;
        return true;
    }

    return false;
}


Vec spatalgs::project(
        const Vec& point,
        const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2)
{
    return plane_p0 + (point - plane_p0).project(plane_p1 - plane_p0, plane_p2 - plane_p0);
}




bool spatalgs::doesRayIntersectPlane(
        const Vec& dir,
        const Vec& pl_p0, const Vec& pl_p1, const Vec& pl_p2)
{
    Vec edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    real_t det = Vec::dot(edges[0], pvec);
    return det > static_cast<real_t>(1e-6) || det < static_cast<real_t>(-1e-6);
}


bool spatalgs::rayIntersectPlane(
    Vec& out_intersectPoint,
    const Vec& origin, const Vec& dir,
    const Vec& pl_p0, const Vec& pl_p1, const Vec& pl_p2)
{
    Vec edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    real_t det = Vec::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    Vec tvec = origin - pl_p0;
    Vec qvec = Vec::cross(tvec, edges[0]);

    real_t t = Vec::dot(edges[1], qvec) / det;
    out_intersectPoint = origin + t * dir;
    return t > static_cast<real_t>(0.0);
}


bool spatalgs::doesRayIntersectTriangle(
    const Vec& origin, const Vec& dir,
    const Vec& tr_p0, const Vec& tr_p1, const Vec& tr_p2)
{
    Vec edges[2]{ tr_p1 - tr_p0, tr_p2 - tr_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    real_t det = Vec::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    real_t inv_det = static_cast<real_t>(1.0) / det;

    Vec tvec = origin - tr_p0;
    real_t u = Vec::dot(tvec, pvec) * inv_det;
    if (u < static_cast<real_t>(0.0) || u > static_cast<real_t>(1.0))
        return false;

    Vec qvec = Vec::cross(tvec, edges[0]);
    real_t v = Vec::dot(dir, qvec) * inv_det;
    if (v < static_cast<real_t>(0.0) || u + v > static_cast<real_t>(1.0))
        return false;

    real_t t = Vec::dot(edges[1], qvec) * inv_det;
    return t > static_cast<real_t>(0.0);
}


Vec spatalgs::linesClosestPoint(const Vec& line0_p0, const Vec& line0_p1, const Vec& line1_p0, const Vec& line1_p1)
{
    Vec u = line0_p1 - line0_p0;
    Vec v = line1_p1 - line1_p0;
    Vec w = line0_p0 - line1_p0;
    real_t a = Vec::dot(u, u); // >= 0
    real_t b = Vec::dot(u, v);
    real_t c = Vec::dot(v, v); // >= 0
    real_t d = Vec::dot(u, w);
    real_t e = Vec::dot(v, w);
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

    return static_cast<real_t>(0.5) * (line0_p0 + sc * u + line1_p0 + tc * v);
}


bool spatalgs::lineIntersectPlane(
    Vec& out_intersectPoint,
    const Vec& line_point, const Vec& line_dir,
    const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2)
{
    Vec edges[2]{ plane_p1 - plane_p0, plane_p2 - plane_p0 };

    Vec pvec = Vec::cross(line_dir, edges[1]);
    real_t det = Vec::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    Vec tvec = line_point - plane_p0;
    Vec qvec = Vec::cross(tvec, edges[0]);

    real_t t = Vec::dot(edges[1], qvec) / det;
    out_intersectPoint = line_point + t * line_dir;
    return true;
}


Vec spatalgs::lineIntersectPlane(
    const Vec& line_point, const Vec& line_dir,
    const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2)
{
    Vec edges[2]{ plane_p1 - plane_p0, plane_p2 - plane_p0 };

    Vec pvec = Vec::cross(line_dir, edges[1]);
    Vec tvec = line_point - plane_p0;
    Vec qvec = Vec::cross(tvec, edges[0]);

    real_t t = Vec::dot(edges[1], qvec) / Vec::dot(edges[0], pvec);
    return line_point + t * line_dir;
}


bool spatalgs::doesSegmentIntersectTriangle(
    const Vec& segm_p0, const Vec& segm_p1,
    const Vec& tr_p0, const Vec& tr_p1, const Vec& tr_p2)
{
    const Vec dir = segm_p1 - segm_p0;

    Vec edges[2]
    { tr_p1 - tr_p0,
        tr_p2 - tr_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    real_t det = Vec::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    real_t inv_det = static_cast<real_t>(1.0) / det;

    Vec tvec = segm_p0 - tr_p0;
    real_t u = Vec::dot(tvec, pvec) * inv_det;
    if (u < static_cast<real_t>(0.0) || u > static_cast<real_t>(1.0))
        return false;

    Vec qvec = Vec::cross(tvec, edges[0]);
    real_t v = Vec::dot(dir, qvec) * inv_det;
    if (v < static_cast<real_t>(0.0) || u + v > static_cast<real_t>(1.0))
        return false;

    real_t t = Vec::dot(edges[1], qvec) * inv_det;
    return t < static_cast<real_t>(1.0) && t > static_cast<real_t>(0.0);
}


bool spatalgs::segmentIntersectPlane(
    Vec& out_intersectPoint,
    const Vec& p0, const Vec& p1,
    const Vec& pl_p0, const Vec& pl_p1, const Vec& pl_p2)
{
    const Vec dir = p1 - p0;

    Vec edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    real_t det = Vec::dot(edges[0], pvec);
    if (det < static_cast<real_t>(1e-6) && det > static_cast<real_t>(-1e-6))
        return false;

    Vec tvec = p0 - pl_p0;
    Vec qvec = Vec::cross(tvec, edges[0]);

    real_t t = Vec::dot(edges[1], qvec) / det;
    out_intersectPoint = p0 + t * dir;
    return t < static_cast<real_t>(1.0) && t > static_cast<real_t>(0.0);
}




bool spatalgs::isPointOnTriangle(
        const Vec& point,
        const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2)
{
    real_t s0 = Vec::cross(trngl_p0 - point,    trngl_p1 - point   ).magnitude();
    real_t s1 = Vec::cross(trngl_p0 - point,    trngl_p2 - point   ).magnitude();
    real_t s2 = Vec::cross(trngl_p1 - point,    trngl_p2 - point   ).magnitude();
    real_t s  = Vec::cross(trngl_p0 - trngl_p2, trngl_p1 - trngl_p2).magnitude();

    real_t expr = s - s0 - s1 - s2;
    return expr > static_cast<real_t>(-1e-6) && expr < static_cast<real_t>(1e-6);
}


bool spatalgs::isPointOnTriangle(
        const Vec& point,
        const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2,
        real_t max_sqrs_sum)
{
    if (computeSqrsSum(point, trngl_p0, trngl_p1, trngl_p2) > max_sqrs_sum)
        return false;

    return isPointOnTriangle(point, trngl_p0, trngl_p1, trngl_p2);
}




real_t spatalgs::linesDistance(
        const Vec& line0_p0, const Vec& line0_p1,
        const Vec& line1_p0, const Vec& line1_p1)
{
    Vec u = line0_p1 - line0_p0;
    Vec v = line1_p1 - line1_p0;
    Vec w = line0_p0 - line1_p0;
    real_t a = Vec::dot(u, u); // >= 0
    real_t b = Vec::dot(u, v);
    real_t c = Vec::dot(v, v); // >= 0
    real_t d = Vec::dot(u, w);
    real_t e = Vec::dot(v, w);
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

    Vec diff_p = w + (sc * u) - (tc * v);

    return diff_p.magnitude();
}


real_t spatalgs::segmentsDistance(
    const Vec& segm0_p0, const Vec& segm0_p1,
    const Vec& segm1_p0, const Vec& segm1_p1)
{
    Vec u = segm0_p1 - segm0_p0;
    Vec v = segm1_p1 - segm1_p0;
    Vec w = segm0_p0 - segm1_p0;
    real_t a = Vec::dot(u, u); // >= 0
    real_t b = Vec::dot(u, v);
    real_t c = Vec::dot(v, v); // >= 0
    real_t d = Vec::dot(u, w);
    real_t e = Vec::dot(v, w);
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

    Vec diff_p = w + (sc * u) - (tc * v);

    return diff_p.magnitude();
}




real_t spatalgs::cpaTime(
    const Vec& start0, const Vec& vel0,
    const Vec& start1, const Vec& vel1)
{
    Vec dv = vel0 - vel1;
    real_t dv2 = Vec::dot(dv, dv);
    if (dv2 < static_cast<real_t>(1e-6))
        return static_cast<real_t>(0.0);

    Vec w0 = start0 - start1;
    return -Vec::dot(w0, dv) / dv2;
}



real_t spatalgs::cpaDistance(
    const Vec& start0, const Vec& vel0,
    const Vec& start1, const Vec& vel1)
{
    real_t cpa_time = cpaTime(start0, vel0, start1, vel1);
    Vec p0 = start0 + cpa_time * vel0;
    Vec p1 = start1 + cpa_time * vel1;
    return (p1 - p0).magnitude();
}




real_t spatalgs::distancePointToLine(const Vec& point, const Vec& line_p0, const Vec& line_p1)
{
    return (project(point, line_p0, line_p1) - point).magnitude();
}


real_t spatalgs::distancePointToSegment(const Vec& point, const Vec& segm_p0, const Vec& segm_p1)
{
    Vec proj = project(point, segm_p0, segm_p1);

    if (INSIDE_RECTANGLE(segm_p0, segm_p1, proj))
    {
        return (proj - point).magnitude();
    }
    else if (real_t sqr_magns[2] { (segm_p0 - point).sqrMagnitude(), (segm_p1 - point).sqrMagnitude() };
             sqr_magns[0] < sqr_magns[1])
    {
        if constexpr (std::is_same<real_t, float>())
        {
            return sqrtf(sqr_magns[0]);
        }
        else
        {
            return sqrt(sqr_magns[0]);
        }
    }
    else
    {
        if constexpr (std::is_same<real_t, float>())
        {
            return sqrtf(sqr_magns[1]);
        }
        else
        {
            return sqrt(sqr_magns[1]);
        }
    }
}


real_t spatalgs::distancePointToTriangleOnPlane(
        const Vec& point, const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2)
{
    Vec closest_points[3];
    closest_points[0] = closestSegmentPointToPoint(point, trngl_p0, trngl_p1);
    closest_points[1] = closestSegmentPointToPoint(point, trngl_p1, trngl_p2);
    closest_points[2] = closestSegmentPointToPoint(point, trngl_p2, trngl_p0);

    real_t sqrs[3];
    sqrs[0] = (closest_points[0] - point).sqrMagnitude();
    sqrs[1] = (closest_points[1] - point).sqrMagnitude();
    sqrs[2] = (closest_points[2] - point).sqrMagnitude();

    return std::min({ sqrs[0], sqrs[1], sqrs[2] });
}


bool spatalgs::doesTriangleIntersectSphere(const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2,
                                           const Vec& center,   real_t radius)
{
    Vec proj = project(center, trngl_p0, trngl_p1, trngl_p2);
    if ((proj - center).sqrMagnitude() > radius * radius)
        return false;

    if (isPointOnTriangle(proj, trngl_p0, trngl_p1, trngl_p2, computeMaxSqrsSum(trngl_p0, trngl_p1, trngl_p2)))
        return true;

    Vec closest = closestTrianglePointToPointOnPlane(proj, trngl_p0, trngl_p1, trngl_p2);
    return (closest - center).sqrMagnitude() <= radius * radius;
}


Vec spatalgs::closestSegmentPointToPoint(const Vec& point, const Vec& segm_p0, const Vec& segm_p1)
{
    Vec proj = project(point, segm_p0, segm_p1);

    if (INSIDE_RECTANGLE(segm_p0, segm_p1, proj))
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


real_t spatalgs::computeMaxSqrsSum(const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2)
{
    real_t sqrs[3];
    sqrs[0] = (trngl_p1 - trngl_p0).sqrMagnitude();
    sqrs[1] = (trngl_p2 - trngl_p1).sqrMagnitude();
    sqrs[2] = (trngl_p0 - trngl_p2).sqrMagnitude();

    int max_inds[2];
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


real_t spatalgs::computeSqrsSum(const Vec& point, const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2)
{
    real_t sqrs[3];
    sqrs[0] = (trngl_p0 - point).sqrMagnitude();
    sqrs[1] = (trngl_p1 - point).sqrMagnitude();
    sqrs[2] = (trngl_p2 - point).sqrMagnitude();
    return sqrs[0] + sqrs[1] + sqrs[2];
}


Vec spatalgs::closestTrianglePointToPointOnPlane(
        const Vec& point, const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2)
{
    Vec closest_points[3];
    closest_points[0] = closestSegmentPointToPoint(point, trngl_p0, trngl_p1);
    closest_points[1] = closestSegmentPointToPoint(point, trngl_p1, trngl_p2);
    closest_points[2] = closestSegmentPointToPoint(point, trngl_p2, trngl_p0);

    real_t sqrs[3];
    sqrs[0] = (closest_points[0] - point).sqrMagnitude();
    sqrs[1] = (closest_points[1] - point).sqrMagnitude();
    sqrs[2] = (closest_points[2] - point).sqrMagnitude();

    int min_i = -1;
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
