// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <algorithm>
#include "helpers/spatial-algs/vec.h"


namespace spatalgs {

Vec  project( const Vec& point,
              const Vec& line_p0, const Vec& line_p1 );

bool project( Vec& out,
              const Vec& point,
              const Vec& segm_p0, const Vec& segm_p1 );

Vec  project( const Vec& point,
                   const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );

bool project( Vec& out,
              const Vec& point,
              const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );


bool doesRayIntersectPlane( const Vec& dir,
                            const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );

bool rayIntersectPlane( Vec& out_intersectPoint,
                        const Vec& origin,   const Vec& dir,
                        const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );


bool doesRayIntersectTriangle( const Vec& origin,   const Vec& dir,
                               const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );

bool rayIntersectTriangle( Vec& out_intersectPoint,
                           const Vec& origin, const Vec& dir,
                           const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );


bool doesLineIntersectPlane( const Vec& line_point, const Vec& line_dir,
                             const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );

bool lineIntersectPlane( Vec& out_intersectPoint,
                         const Vec& line_point, const Vec& line_dir,
                         const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );

Vec  lineIntersectPlane( const Vec& line_point, const Vec& line_dir,
                               const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );


bool doesSegmentIntersectTriangle( const Vec& segm_p0, const Vec& segm_p1,
                                   const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );

bool segmentIntersectTriangle( Vec& out_intersectPoint,
                               const Vec& segm_p0,  const Vec& segm_p1,
                               const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );


bool doesSegmentIntersectPlane( const Vec& segm_p0,  const Vec& segm_p1,
                                const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );

bool segmentIntersectPlane( Vec& out_intersectPoint,
                            const Vec& segm_p0,  const Vec& segm_p1,
                            const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );

Vec segmentIntersectPlane( const Vec& segm_p0,  const Vec& segm_p1,
                                 const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );


bool doesTriangleIntersectSphere( const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2,
                                  const Vec& center,   real_t radius );


real_t computeSqrsSum( const Vec& point,
                       const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );

real_t computeMaxSqrsSum( const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );

bool isPointOnPlane(    const Vec& point,
                        const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );

bool isPointOnTriangle( const Vec& point,
                        const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );

bool isPointOnTriangle( const Vec& point,
                        const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2,
                        real_t max_sqrs_sum );


Vec closestSegmentPointToPoint( const Vec& point,
                                const Vec& segm_p0, const Vec& segm_p1 );

Vec closestTrianglePointToPointOnPlane( const Vec& point,
                                        const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );

real_t distancePointToLine(     const Vec& point,
                              const Vec& line_p0,  const Vec& line_p1 );

real_t distancePointToSegment( const Vec& point,
                               const Vec& segm_p0,  const Vec& segm_p1 );

real_t distancePointToPlane( const Vec& point,
                             const Vec& plane_p0, const Vec& plane_p1, const Vec& plane_p2 );

real_t distancePointToTriangle( const Vec& point,
                                const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );

real_t distancePointToTriangleOnPlane( const Vec& point,
                                       const Vec& trngl_p0, const Vec& trngl_p1, const Vec& trngl_p2 );


Vec linesClosestPoint( const Vec& line0_p0, const Vec& line0_p1,
                       const Vec& line1_p0, const Vec& line1_p1 );

real_t linesDistance( const Vec& line0_p0, const Vec& line0_p1,
                      const Vec& line1_p0, const Vec& line1_p1 );

real_t segmentsDistance( const Vec& segm0_p0, const Vec& segm0_p1,
                         const Vec& segm1_p0, const Vec& segm1_p1 );


// CPA - Closest Vec of Approach.
real_t cpaTime( const Vec& start0, const Vec& vel0,
                const Vec& start1, const Vec& vel1 );

real_t cpaDistance( const Vec& start0, const Vec& vel0,
                    const Vec& start1, const Vec& vel1 );

} // namespace spatalgs
