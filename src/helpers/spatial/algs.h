// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <algorithm>
#include "vec.h"

// TODO: experiment with pass by value instead of reference and make benchmark
namespace spt::algs {

vec3 project(const vec3& point,
             const vec3& line_p0, const vec3& line_p1);

// TODO: return std::optional<vec3>
bool project(vec3& out,
             const vec3& point,
             const vec3& segm_p0, const vec3& segm_p1);

vec3 project(const vec3& point,
             const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);

// TODO: return std::optional<vec3>
bool project(vec3& out,
             const vec3& point,
             const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);


bool doesRayIntersectPlane(const vec3& dir,
                           const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);

// TODO: return std::optional<vec3>
bool rayIntersectPlane(vec3& out_intersectPoint,
                       const vec3& origin, const vec3& dir,
                       const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);


bool doesRayIntersectTriangle(const vec3& origin, const vec3& dir,
                              const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);

// TODO: return std::optional<vec3>
bool rayIntersectTriangle(vec3& out_intersectPoint,
                          const vec3& origin, const vec3& dir,
                          const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);


bool doesLineIntersectPlane(const vec3& line_point, const vec3& line_dir,
                            const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);

// TODO: return std::optional<vec3>
bool lineIntersectPlane(vec3& out_intersectPoint,
                        const vec3& line_point, const vec3& line_dir,
                        const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);

vec3 lineIntersectPlane(const vec3& line_point, const vec3& line_dir,
                        const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);


bool doesSegmentIntersectTriangle(const vec3& segm_p0, const vec3& segm_p1,
                                  const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);

// TODO: return std::optional<vec3>
bool segmentIntersectTriangle(vec3& out_intersectPoint,
                              const vec3& segm_p0, const vec3& segm_p1,
                              const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);


bool doesSegmentIntersectPlane(const vec3& segm_p0, const vec3& segm_p1,
                               const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);

// TODO: return std::optional<vec3>
bool segmentIntersectPlane(vec3& out_intersectPoint,
                           const vec3& segm_p0, const vec3& segm_p1,
                           const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);

vec3 segmentIntersectPlane(const vec3& segm_p0, const vec3& segm_p1,
                           const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);


bool doesTriangleIntersectSphere(const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2,
                                 const vec3& center, real_t radius);

real_t computeSqrsSum(const vec3& point,
                      const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);

real_t computeMaxSqrsSum(const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);

bool isPointOnPlane(const vec3& point,
                    const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);

bool isPointOnTriangle(const vec3& point,
                       const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);

bool isPointOnTriangle(const vec3& point,
                       const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2,
                       real_t max_sqrs_sum);

bool isPointInTetrahedron(const vec3& point,
                          const vec3& tetr_p0, const vec3& tetr_p1, const vec3& tetr_p2, const vec3& tetr_p3);

bool isPointInCylinder(const vec3& point,
                       const vec3& cyl_p0, const vec3& cyl_p1, real_t radius);

vec3 closestSegmentPointToPoint(const vec3& point,
                                const vec3& segm_p0, const vec3& segm_p1);

vec3 closestTrianglePointToPointOnPlane(const vec3& point,
                                        const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);

real_t distancePointToLine(const vec3& point,
                           const vec3& line_p0, const vec3& line_p1);

real_t distancePointToSegment(const vec3& point,
                              const vec3& segm_p0, const vec3& segm_p1);

real_t distancePointToPlane(const vec3& point,
                            const vec3& plane_p0, const vec3& plane_p1, const vec3& plane_p2);

real_t distancePointToTriangle(const vec3& point,
                               const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);

real_t distancePointToTriangleOnPlane(const vec3& point,
                                      const vec3& trngl_p0, const vec3& trngl_p1, const vec3& trngl_p2);

vec3 linesClosestPoint(const vec3& line0_p0, const vec3& line0_p1,
                       const vec3& line1_p0, const vec3& line1_p1);

real_t linesDistance(const vec3& line0_p0, const vec3& line0_p1,
                     const vec3& line1_p0, const vec3& line1_p1);

real_t segmentsDistance(const vec3& segm0_p0, const vec3& segm0_p1,
                        const vec3& segm1_p0, const vec3& segm1_p1);

// CPA - Closest vec3 of Approach.
real_t cpaTime(const vec3& start0, const vec3& vel0,
               const vec3& start1, const vec3& vel1);

real_t cpaDistance(const vec3& start0, const vec3& vel0,
                   const vec3& start1, const vec3& vel1);

} // namespace spt::algs
