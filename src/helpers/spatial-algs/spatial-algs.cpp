// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#include "helpers/spatial-algs/spatial-algs.h"

using tva::Vec;
using tva::Point;
namespace spatalgs = tva::spatalgs;




#define DET(a, b, c, d) \
        ((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0.coors[0], corner1.coors[0], point.coors[0]) && \
         BETWEEN(corner0.coors[1], corner1.coors[1], point.coors[1]) && \
         BETWEEN(corner0.coors[2], corner1.coors[2], point.coors[2]))





Point spatalgs::project(
        const Point& point,
        const Point& line_p0, const Point& line_p1)
{
    return line_p0 + (point - line_p0).project(line_p1 - line_p0);
}


bool spatalgs::project(
    Point& out,
    const Point& point,
    const Point& segm_p0, const Point& segm_p1)
{
    Vec res = segm_p0 + (point - segm_p0).project(segm_p1 - segm_p0);

    const double EPS = (segm_p1 - segm_p0).magnitude() * 1e-6;
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


Point spatalgs::project(
        const Point& point,
        const Point& plane_p0, const Point& plane_p1, const Point& plane_p2)
{
    return plane_p0 + (point - plane_p0).project(plane_p1 - plane_p0, plane_p2 - plane_p0);
}




bool spatalgs::doesRayIntersectPlane(
        const Vec& dir,
        const Point& pl_p0, const Point& pl_p1, const Point& pl_p2)
{
    Vec edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    double det = Vec::dot(edges[0], pvec);
    return det > 1e-6 || det < -1e-6;
}


bool spatalgs::rayIntersectPlane(
    Point& out_intersectPoint,
    const Point& origin, const Vec& dir,
    const Point& pl_p0, const Point& pl_p1, const Point& pl_p2)
{
    Vec edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    double det = Vec::dot(edges[0], pvec);
    if (det < 1e-6 && det > -1e-6)
        return false;

    Vec tvec = origin - pl_p0;
    Vec qvec = Vec::cross(tvec, edges[0]);

    double t = Vec::dot(edges[1], qvec) / det;
    out_intersectPoint = origin + t * dir;
    return t > 0.0;
}


bool spatalgs::doesRayIntersectTriangle(
    const Point& origin, const Vec& dir,
    const Point& tr_p0, const Point& tr_p1, const Point& tr_p2)
{
    Vec edges[2]{ tr_p1 - tr_p0, tr_p2 - tr_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    double det = Vec::dot(edges[0], pvec);
    if (det < 1e-6 && det > -1e-6)
        return false;

    double inv_det = 1.0 / det;

    Vec tvec = origin - tr_p0;
    double u = Vec::dot(tvec, pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    Vec qvec = Vec::cross(tvec, edges[0]);
    double v = Vec::dot(dir, qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    double t = Vec::dot(edges[1], qvec) * inv_det;
    return t > 0.0;
}


tva::Point spatalgs::linesClosestPoint(const tva::Point& line0_p0, const tva::Point& line0_p1, const tva::Point& line1_p0, const tva::Point& line1_p1)
{
    Vec u = line0_p1 - line0_p0;
    Vec v = line1_p1 - line1_p0;
    Vec w = line0_p0 - line1_p0;
    double a = Vec::dot(u, u); // >= 0
    double b = Vec::dot(u, v);
    double c = Vec::dot(v, v); // >= 0
    double d = Vec::dot(u, w);
    double e = Vec::dot(v, w);
    double determ = DET(a, b, b, c); // >= 0
    double sc, tc;

    if (determ < 1e-6)
    {
        sc = 0.0;
        tc = b > c ? d / b : e / c;
    }
    else
    {
        double inv_determ = 1.0 / determ;
        sc = DET(b, c, d, e) * inv_determ;
        tc = DET(a, b, d, e) * inv_determ;
    }

    return 0.5 * (line0_p0 + sc * u + line1_p0 + tc * v);
}


bool spatalgs::lineIntersectPlane(
    Point& out_intersectPoint,
    const Point& line_point, const Point& line_dir,
    const Point& plane_p0, const Point& plane_p1, const Point& plane_p2)
{
    Vec edges[2]{ plane_p1 - plane_p0, plane_p2 - plane_p0 };

    Vec pvec = Vec::cross(line_dir, edges[1]);
    double det = Vec::dot(edges[0], pvec);
    if (det < 1e-6 && det > -1e-6)
        return false;

    Vec tvec = line_point - plane_p0;
    Vec qvec = Vec::cross(tvec, edges[0]);

    double t = Vec::dot(edges[1], qvec) / det;
    out_intersectPoint = line_point + t * line_dir;
    return true;
}


Point spatalgs::lineIntersectPlane(
    const Point& line_point, const Point& line_dir,
    const Point& plane_p0, const Point& plane_p1, const Point& plane_p2)
{
    Vec edges[2]{ plane_p1 - plane_p0, plane_p2 - plane_p0 };

    Vec pvec = Vec::cross(line_dir, edges[1]);
    Vec tvec = line_point - plane_p0;
    Vec qvec = Vec::cross(tvec, edges[0]);

    double t = Vec::dot(edges[1], qvec) / Vec::dot(edges[0], pvec);
    return line_point + t * line_dir;
}


bool spatalgs::doesSegmentIntersectTriangle(
    const Point& segm_p0, const Point& segm_p1,
    const Point& tr_p0, const Point& tr_p1, const Point& tr_p2)
{
    const Vec dir = segm_p1 - segm_p0;

    Vec edges[2]
    { tr_p1 - tr_p0,
        tr_p2 - tr_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    double det = Vec::dot(edges[0], pvec);
    if (det < 1e-6 && det > -1e-6)
        return false;

    double inv_det = 1.0 / det;

    Vec tvec = segm_p0 - tr_p0;
    double u = Vec::dot(tvec, pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    Vec qvec = Vec::cross(tvec, edges[0]);
    double v = Vec::dot(dir, qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    double t = Vec::dot(edges[1], qvec) * inv_det;
    return t < 1.0 && t > 0.0;
}


bool spatalgs::segmentIntersectPlane(
    Point& out_intersectPoint,
    const Point& p0, const Point& p1,
    const Point& pl_p0, const Point& pl_p1, const Point& pl_p2)
{
    const Vec dir = p1 - p0;

    Vec edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

    Vec pvec = Vec::cross(dir, edges[1]);
    double det = Vec::dot(edges[0], pvec);
    if (det < 1e-6 && det > -1e-6)
        return false;

    Vec tvec = p0 - pl_p0;
    Vec qvec = Vec::cross(tvec, edges[0]);

    double t = Vec::dot(edges[1], qvec) / det;
    out_intersectPoint = p0 + t * dir;
    return t < 1.0 && t > 0.0;
}




bool spatalgs::isPointOnTriangle(
        const Point& point,
        const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2)
{
    double s0 = Vec::cross(trngl_p0 - point,    trngl_p1 - point   ).magnitude();
    double s1 = Vec::cross(trngl_p0 - point,    trngl_p2 - point   ).magnitude();
    double s2 = Vec::cross(trngl_p1 - point,    trngl_p2 - point   ).magnitude();
    double s  = Vec::cross(trngl_p0 - trngl_p2, trngl_p1 - trngl_p2).magnitude();

    double expr = s - s0 - s1 - s2;
    return expr > -1e-6 && expr < 1e-6;
}


bool spatalgs::isPointOnTriangle(
        const Point& point,
        const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2,
        double max_sqrs_sum)
{
    if (computeSqrsSum(point, trngl_p0, trngl_p1, trngl_p2) > max_sqrs_sum)
        return false;

    return isPointOnTriangle(point, trngl_p0, trngl_p1, trngl_p2);
}




double spatalgs::linesDistance(
        const Point& line0_p0, const Point& line0_p1,
        const Point& line1_p0, const Point& line1_p1)
{
    Vec u = line0_p1 - line0_p0;
    Vec v = line1_p1 - line1_p0;
    Vec w = line0_p0 - line1_p0;
    double a = Vec::dot(u, u); // >= 0
    double b = Vec::dot(u, v);
    double c = Vec::dot(v, v); // >= 0
    double d = Vec::dot(u, w);
    double e = Vec::dot(v, w);
    double determ = DET(a, b, b, c); // >= 0
    double sc, tc;

    if (determ < 1e-6)
    {
        sc = 0.0;
        tc = b > c ? d / b : e / c;
    }
    else
    {
        double inv_determ = 1.0 / determ;
        sc = DET(b, c, d, e) * inv_determ;
        tc = DET(a, b, d, e) * inv_determ;
    }

    Vec diff_p = w + (sc * u) - (tc * v);

    return diff_p.magnitude();
}


double spatalgs::segmentsDistance(
    const Point& segm0_p0, const Point& segm0_p1,
    const Point& segm1_p0, const Point& segm1_p1)
{
    Vec u = segm0_p1 - segm0_p0;
    Vec v = segm1_p1 - segm1_p0;
    Vec w = segm0_p0 - segm1_p0;
    double a = Vec::dot(u, u); // >= 0
    double b = Vec::dot(u, v);
    double c = Vec::dot(v, v); // >= 0
    double d = Vec::dot(u, w);
    double e = Vec::dot(v, w);
    double determ = DET(a, b, b, c); // >= 0
    double sc, sn, sd = determ;
    double tc, tn, td = determ;

    if (determ < 1e-6)
    {
        sn = 0.0;
        sd = 1.0;
        tn = e;
        td = c;
    }
    else
    {
        sn = DET(b, c, d, e);
        tn = DET(a, b, d, e);

        if (sn < 0.0)
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

    if (tn < 0.0)
    {
        tn = 0.0;

        if (-d < 0.0)
        {
            sn = 0.0;
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

        if (b - d < 0.0)
        {
            sn = 0;
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

    sc = std::abs(sn) < 1e-6 ? 0.0 : sn / sd;
    tc = std::abs(tn) < 1e-6 ? 0.0 : tn / td;

    Vec diff_p = w + (sc * u) - (tc * v);

    return diff_p.magnitude();
}




double spatalgs::cpaTime(
    const Point& start0, const Vec& vel0,
    const Point& start1, const Vec& vel1)
{
    Vec dv = vel0 - vel1;
    double dv2 = Vec::dot(dv, dv);
    if (dv2 < 1e-6)
        return 0.0;

    Vec w0 = start0 - start1;
    return -Vec::dot(w0, dv) / dv2;
}


double spatalgs::cpaDistance(
    const Point& start0, const Vec& vel0,
    const Point& start1, const Vec& vel1)
{
    double cpa_time = cpaTime(start0, vel0, start1, vel1);
    Vec p0 = start0 + cpa_time * vel0;
    Vec p1 = start1 + cpa_time * vel1;
    return (p1 - p0).magnitude();
}




double spatalgs::distancePointToLine(const Point& point, const Point& line_p0, const Point& line_p1)
{
    return (project(point, line_p0, line_p1) - point).magnitude();
}


double spatalgs::distancePointToSegment(const Point& point, const Point& segm_p0, const Point& segm_p1)
{
    Point proj = project(point, segm_p0, segm_p1);

    // Try to remove later.
    const double EPS = (segm_p1 - segm_p0).magnitude() * 1e-6;
    if (INSIDE_RECTANGLE(segm_p0, segm_p1, proj))
    {
        return (proj - point).magnitude();
    }
    else if (double sqr_magns[2] { (segm_p0 - point).sqrMagnitude(), (segm_p1 - point).sqrMagnitude() };
             sqr_magns[0] < sqr_magns[1])
    {
        return sqrt(sqr_magns[0]);
    }
    else
    {
        return sqrt(sqr_magns[1]);
    }
}


double spatalgs::distancePointToTriangleOnPlane(
        const Point& point, const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2)
{
    Point closest_points[3];
    closest_points[0] = closestSegmentPointToPoint(point, trngl_p0, trngl_p1);
    closest_points[1] = closestSegmentPointToPoint(point, trngl_p1, trngl_p2);
    closest_points[2] = closestSegmentPointToPoint(point, trngl_p2, trngl_p0);

    double sqrs[3];
    sqrs[0] = (closest_points[0] - point).sqrMagnitude();
    sqrs[1] = (closest_points[1] - point).sqrMagnitude();
    sqrs[2] = (closest_points[2] - point).sqrMagnitude();

    return std::min({ sqrs[0], sqrs[1], sqrs[2] });
}


bool spatalgs::doesTriangleIntersectSphere(const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2,
                                           const Point& center,   double radius)
{
    Point proj = project(center, trngl_p0, trngl_p1, trngl_p2);
    if ((proj - center).sqrMagnitude() > radius * radius)
        return false;

    if (isPointOnTriangle(proj, trngl_p0, trngl_p1, trngl_p2, computeMaxSqrsSum(trngl_p0, trngl_p1, trngl_p2)))
        return true;

    Point closest = closestTrianglePointToPointOnPlane(proj, trngl_p0, trngl_p1, trngl_p2);
    return (closest - center).sqrMagnitude() <= radius * radius;
}


Point spatalgs::closestSegmentPointToPoint(const Point& point, const Point& segm_p0, const Point& segm_p1)
{
    Point proj = project(point, segm_p0, segm_p1);

    // Try to remove later.
    const double EPS = (segm_p1 - segm_p0).magnitude() * 1e-6;
    if (INSIDE_RECTANGLE(segm_p0, segm_p1, proj))
    {
        return proj;
    }
    else if (double sqr_magns[2] { (segm_p0 - point).sqrMagnitude(), (segm_p1 - point).sqrMagnitude() };
             sqr_magns[0] < sqr_magns[1])
    {
        return segm_p0;
    }
    else
    {
        return segm_p1;
    }
}


double spatalgs::computeMaxSqrsSum(const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2)
{
    double sqrs[3];
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


double spatalgs::computeSqrsSum(const Point& point, const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2)
{
    double sqrs[3];
    sqrs[0] = (trngl_p0 - point).sqrMagnitude();
    sqrs[1] = (trngl_p1 - point).sqrMagnitude();
    sqrs[2] = (trngl_p2 - point).sqrMagnitude();
    return sqrs[0] + sqrs[1] + sqrs[2];
}


Point spatalgs::closestTrianglePointToPointOnPlane(
        const Point& point, const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2)
{
    Point closest_points[3];
    closest_points[0] = closestSegmentPointToPoint(point, trngl_p0, trngl_p1);
    closest_points[1] = closestSegmentPointToPoint(point, trngl_p1, trngl_p2);
    closest_points[2] = closestSegmentPointToPoint(point, trngl_p2, trngl_p0);

    double sqrs[3];
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
