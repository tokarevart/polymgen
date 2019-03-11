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

#define PI         3.141592653589793
#define PI_DIV_180 0.017453292519943295




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
    double a = Vec::dot(u, u); // Always >= 0
    double b = Vec::dot(u, v);
    double c = Vec::dot(v, v); // Always >= 0
    double d = Vec::dot(u, w);
    double e = Vec::dot(v, w);
    double D = DET(a, b, b, c); // Always >= 0
    double sc, tc;

    // Compute the line parameters of the two closest points
    if (D < 1e-6) // The lines are almost parallel
    {
        sc = 0.0;
        tc = (b > c ? d / b : e / c); // Use the largest denominator
    }
    else
    {
        double inv_D = 1.0 / D;
        sc = DET(b, c, d, e) * inv_D;
        tc = DET(a, b, d, e) * inv_D;
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


const Point spatalgs::lineIntersectPlane(
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
    double s0 = Vec::cross(trngl_p0 - point, trngl_p1 - point).magnitude();
    double s1 = Vec::cross(trngl_p0 - point, trngl_p2 - point).magnitude();
    double s2 = Vec::cross(trngl_p1 - point, trngl_p2 - point).magnitude();
    double s = Vec::cross(trngl_p0 - trngl_p2, trngl_p1 - trngl_p2).magnitude();

    return (s - s0 - s1 - s2) > -1e-6 && (s - s0 - s1 - s2) < 1e-6;
}




double spatalgs::linesDistance(
    const Point& line0_p0, const Point& line0_p1,
    const Point& line1_p0, const Point& line1_p1)
{
    Vec u = line0_p1 - line0_p0;
    Vec v = line1_p1 - line1_p0;
    Vec w = line0_p0 - line1_p0;
    double a = Vec::dot(u, u); // Always >= 0
    double b = Vec::dot(u, v);
    double c = Vec::dot(v, v); // Always >= 0
    double d = Vec::dot(u, w);
    double e = Vec::dot(v, w);
    double D = DET(a, b, b, c); // Always >= 0
    double sc, tc;

    // Compute the line parameters of the two closest points
    if (D < 1e-6) // The lines are almost parallel
    {
        sc = 0.0;
        tc = (b > c ? d / b : e / c); // Use the largest denominator
    }
    else
    {
        double inv_D = 1.0 / D;
        sc = DET(b, c, d, e) * inv_D;
        tc = DET(a, b, d, e) * inv_D;
    }

    // Get the difference of the two closest points
    Vec diff_p = w + (sc * u) - (tc * v);  // =  L1(sc) - L2(tc)

    return diff_p.magnitude();   // Return the closest distance
}


double spatalgs::segmentsDistance(
    const Point& segm0_p0, const Point& segm0_p1,
    const Point& segm1_p0, const Point& segm1_p1)
{
    Vec u = segm0_p1 - segm0_p0;
    Vec v = segm1_p1 - segm1_p0;
    Vec w = segm0_p0 - segm1_p0;
    double a = Vec::dot(u, u); // Always >= 0
    double b = Vec::dot(u, v);
    double c = Vec::dot(v, v); // Always >= 0
    double d = Vec::dot(u, w);
    double e = Vec::dot(v, w);
    double D = DET(a, b, b, c); // Always >= 0
    double sc, sN, sD = D; // sc = sN / sD, default sD = D >= 0
    double tc, tN, tD = D; // tc = tN / tD, default tD = D >= 0

                           // Compute the line parameters of the two closest points
    if (D < 1e-6) // The lines are almost parallel
    {
        sN = 0.0; // Force using point P0 on segment S1
        sD = 1.0; // to prevent possible division by 0.0 later
        tN = e;
        tD = c;
    }
    else // Get the closest points on the infinite lines
    {
        sN = DET(b, c, d, e);
        tN = DET(a, b, d, e);

        if (sN < 0.0) // sc < 0 => the s=0 shell::Edge is visible
        {
            sN = 0.0;
            tN = e;
            tD = c;
        }
        else if (sN > sD) // sc > 1  => the s=1 shell::Edge is visible
        {
            sN = sD;
            tN = e + b;
            tD = c;
        }
    }

    if (tN < 0.0) // tc < 0 => the t=0 shell::Edge is visible
    {
        tN = 0.0;

        // Recompute sc for this shell::Edge
        if (-d < 0.0)
        {
            sN = 0.0;
        }
        else if (-d > a)
        {
            sN = sD;
        }
        else
        {
            sN = -d;
            sD = a;
        }
    }
    else if (tN > tD) // tc > 1  => the t=1 shell::Edge is visible
    {
        tN = tD;

        // Recompute sc for this shell::Edge
        if (b - d < 0.0)
        {
            sN = 0;
        }
        else if (b - d > a)
        {
            sN = sD;
        }
        else
        {
            sN = b - d;
            sD = a;
        }
    }

    // Finally do the division to get sc and tc
    sc = (std::abs(sN) < 1e-6 ? 0.0 : sN / sD);
    tc = (std::abs(tN) < 1e-6 ? 0.0 : tN / tD);

    // Get the difference of the two closest points
    Vec diff_p = w + (sc * u) - (tc * v); // =  S1(sc) - S2(tc)

    return diff_p.magnitude(); // Return the closest distance
}




double spatalgs::cpaTime(
    const Point& start0, const Vec& vel0,
    const Point& start1, const Vec& vel1)
{
    Vec dv = vel0 - vel1;
    double dv2 = Vec::dot(dv, dv);
    if (dv2 < 1e-6) // The tracks are almost parallel.
        return 0.0;

    Vec w0 = start0 - start1;
    return -Vec::dot(w0, dv) / dv2;
}


double spatalgs::cpaDistance(
    const Point& start0, const Vec& vel0,
    const Point& start1, const Vec& vel1)
{
    double cpa_time = cpaTime(start0, vel0, start1, vel1);
    Vec p0 = start0 + (cpa_time * vel0);
    Vec p1 = start1 + (cpa_time * vel1);
    return (p1 - p0).magnitude();
}
