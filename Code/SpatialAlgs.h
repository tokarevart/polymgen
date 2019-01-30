#pragma once
#include <cmath>
#include <algorithm>
#include "Vec3.h"


namespace spatialalgs
{
    Point3 project(
        const Point3& point,
        const Point3& line_p0, const Point3& line_p1);
    const bool project(
        Point3& out,
        const Point3& point, 
        const Point3& segm_p0, const Point3& segm_p1);

    const bool isRayIntersectPlane(
        const Point3& origin, const Vec3& dir, 
        const Point3& plane_p0, const Point3& plane_p1, const Point3& plane_p2);
    const bool isRayIntersectTriangle(
        const Point3& origin, const Vec3& dir, 
        const Point3& trngl_p0, const Point3& trngl_p1, const Point3& trngl_p2);
    const bool rayIntersectPlane(
        Point3& out_intersectPoint, 
        const Point3& origin, const Vec3& dir, 
        const Point3& plane_p0, const Point3& plane_p1, const Point3& plane_p2);
    const bool rayIntersectTriangle(
        Point3& out_intersectPoint, 
        const Point3& origin, const Vec3& dir, 
        const Point3& trngl_p0, const Point3& trngl_p1, const Point3& trngl_p2);
    const bool isLineIntersectPlane(
        const Point3& origin, const Vec3& dir, 
        const Point3& plane_p0, const Point3& plane_p1, const Point3& plane_p2);
    const bool lineIntersectPlane(
        Point3& out_intersectPoint,
        const Point3& line_point, const Point3& line_dir,
        const Point3& plane_p0, const Point3& plane_p1, const Point3& plane_p2);
    const Point3 lineIntersectPlane(
        const Point3& line_point, const Point3& line_dir,
        const Point3& plane_p0, const Point3& plane_p1, const Point3& plane_p2);
    const bool isSegmentIntersectTriangle(
        const Point3& segm_p0, const Point3& segm_p1, 
        const Point3& trngl_p0, const Point3& trngl_p1, const Point3& trngl_p2);
    const bool isSegmentIntersectPlane(
        const Point3& segm_p0, const Point3& segm_p1, 
        const Point3& plane_p0, const Point3& plane_p1, const Point3& plane_p2);
    const bool segmentIntersectTriangle(
        Point3& out_intersectPoint,
        const Point3& segm_p0, const Point3& segm_p1, 
        const Point3& trngl_p0, const Point3& trngl_p1, const Point3& trngl_p2);
    const bool segmentIntersectPlane(
        Point3& out_intersectPoint,
        const Point3& segm_p0, const Point3& segm_p1, 
        const Point3& plane_p0, const Point3& plane_p1, const Point3& plane_p2);

    const double linesDistance(
        const Point3& line0_p0, const Point3& line0_p1, 
        const Point3& line1_p0, const Point3& line1_p1);
    const double segmentsDistance(
        const Point3& segm0_p0, const Point3& segm0_p1,
        const Point3& segm1_p0, const Point3& segm1_p1);

    // CPA - Closest Point of Approach.
    const double cpaTime(
        const Point3& start0, const Vec3& vel0, 
        const Point3& start1, const Vec3& vel1);
    const double cpaDistance(
        const Point3& start0, const Vec3& vel0, 
        const Point3& start1, const Vec3& vel1);
}