#pragma once
#include <cmath>
#include <algorithm>
#include "Vec3.h"


namespace tva
{
    namespace spatialalgs
    {
        tva::Point3 project(
            const tva::Point3& point,
            const tva::Point3& line_p0, const tva::Point3& line_p1);
        bool project(
            tva::Point3& out,
            const tva::Point3& point,
            const tva::Point3& segm_p0, const tva::Point3& segm_p1);

        bool isRayIntersectPlane(
            const tva::Vec3& dir,
            const tva::Point3& plane_p0, const tva::Point3& plane_p1, const tva::Point3& plane_p2);
        bool isRayIntersectTriangle(
            const tva::Point3& origin, const tva::Vec3& dir,
            const tva::Point3& trngl_p0, const tva::Point3& trngl_p1, const tva::Point3& trngl_p2);
        bool rayIntersectPlane(
            tva::Point3& out_intersectPoint,
            const tva::Point3& origin, const tva::Vec3& dir,
            const tva::Point3& plane_p0, const tva::Point3& plane_p1, const tva::Point3& plane_p2);
        bool rayIntersectTriangle(
            tva::Point3& out_intersectPoint,
            const tva::Point3& origin, const tva::Vec3& dir,
            const tva::Point3& trngl_p0, const tva::Point3& trngl_p1, const tva::Point3& trngl_p2);
        bool isLineIntersectPlane(
            const tva::Point3& origin, const tva::Vec3& dir,
            const tva::Point3& plane_p0, const tva::Point3& plane_p1, const tva::Point3& plane_p2);
        bool lineIntersectPlane(
            tva::Point3& out_intersectPoint,
            const tva::Point3& line_point, const tva::Point3& line_dir,
            const tva::Point3& plane_p0, const tva::Point3& plane_p1, const tva::Point3& plane_p2);
        const tva::Point3 lineIntersectPlane(
            const tva::Point3& line_point, const tva::Point3& line_dir,
            const tva::Point3& plane_p0, const tva::Point3& plane_p1, const tva::Point3& plane_p2);
        bool isSegmentIntersectTriangle(
            const tva::Point3& segm_p0, const tva::Point3& segm_p1,
            const tva::Point3& trngl_p0, const tva::Point3& trngl_p1, const tva::Point3& trngl_p2);
        bool isSegmentIntersectPlane(
            const tva::Point3& segm_p0, const tva::Point3& segm_p1,
            const tva::Point3& plane_p0, const tva::Point3& plane_p1, const tva::Point3& plane_p2);
        bool segmentIntersectTriangle(
            tva::Point3& out_intersectPoint,
            const tva::Point3& segm_p0, const tva::Point3& segm_p1,
            const tva::Point3& trngl_p0, const tva::Point3& trngl_p1, const tva::Point3& trngl_p2);
        bool segmentIntersectPlane(
            tva::Point3& out_intersectPoint,
            const tva::Point3& segm_p0, const tva::Point3& segm_p1,
            const tva::Point3& plane_p0, const tva::Point3& plane_p1, const tva::Point3& plane_p2);

        bool isPointOnTriangle(
            const tva::Point3& point,
            const tva::Point3& trngl_p0, const tva::Point3& trngl_p1, const tva::Point3& trngl_p2);
        bool isPointOnPlane(
            const tva::Point3& point,
            const tva::Point3& plane_p0, const tva::Point3& plane_p1, const tva::Point3& plane_p2);

        double linesDistance(
            const tva::Point3& line0_p0, const tva::Point3& line0_p1,
            const tva::Point3& line1_p0, const tva::Point3& line1_p1);
        double segmentsDistance(
            const tva::Point3& segm0_p0, const tva::Point3& segm0_p1,
            const tva::Point3& segm1_p0, const tva::Point3& segm1_p1);

        // CPA - Closest Point of Approach.
        double cpaTime(
            const tva::Point3& start0, const tva::Vec3& vel0,
            const tva::Point3& start1, const tva::Vec3& vel1);
        double cpaDistance(
            const tva::Point3& start0, const tva::Vec3& vel0,
            const tva::Point3& start1, const tva::Vec3& vel1);
    }
}
