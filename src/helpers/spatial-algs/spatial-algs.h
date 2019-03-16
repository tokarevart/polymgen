#pragma once
#include <cmath>
#include <algorithm>
#include "helpers/spatial-algs/vec.h"


namespace tva {
namespace spatalgs {

Point project( const Point& point,
               const Point& line_p0, const Point& line_p1 );

bool  project( Point& out,
               const Point& point,
               const Point& segm_p0, const Point& segm_p1 );

Point project( const Point& point,
               const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );

bool  project( Point& out,
               const Point& point,
               const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );


bool doesRayIntersectPlane( const Vec& dir,
                            const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );

bool rayIntersectPlane( Point& out_intersectPoint,
                        const Point& origin, const Vec& dir,
                        const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );


bool doesRayIntersectTriangle( const Point& origin, const Vec& dir,
                               const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );

bool rayIntersectTriangle( Point& out_intersectPoint,
                           const Point& origin, const Vec& dir,
                           const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );


bool doesLineIntersectPlane( const Point& line_point, const Point& line_dir,
                             const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );

bool lineIntersectPlane( Point& out_intersectPoint,
                         const Point& line_point, const Point& line_dir,
                         const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );

Point lineIntersectPlane( const Point& line_point, const Point& line_dir,
                          const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );


bool doesSegmentIntersectTriangle( const Point& segm_p0, const Point& segm_p1,
                                   const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );

bool segmentIntersectTriangle( Point& out_intersectPoint,
                               const Point& segm_p0, const Point& segm_p1,
                               const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );


bool doesSegmentIntersectPlane( const Point& segm_p0, const Point& segm_p1,
                                const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );

bool segmentIntersectPlane( Point& out_intersectPoint,
                            const Point& segm_p0, const Point& segm_p1,
                            const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );

Point segmentIntersectPlane( const Point& segm_p0, const Point& segm_p1,
                             const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );


bool doesTriangleIntersectSphere( const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2,
                                  const Point& center,   double radius );


double computeSqrsSum( const Point& point,
                       const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );

double computeMaxSqrsSum( const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );

bool isPointOnPlane(    const Point& point,
                        const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );

bool isPointOnTriangle( const Point& point,
                        const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );

bool isPointOnTriangle( const Point& point,
                        const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2,
                        double max_sqrs_sum );


Point closestSegmentPointToPoint( const Point& point,
                                  const Point& segm_p0, const Point& segm_p1 );

Point closestTrianglePointToPointOnPlane( const Point& point,
                                          const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );

double distancePointToLine(     const Point& point,
                                const Point& line_p0, const Point& line_p1 );

double distancePointToSegment(  const Point& point,
                                const Point& segm_p0, const Point& segm_p1 );

double distancePointToPlane(    const Point& point,
                                const Point& plane_p0, const Point& plane_p1, const Point& plane_p2 );

double distancePointToTriangle( const Point& point,
                                const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );

double distancePointToTriangleOnPlane( const Point& point,
                                       const Point& trngl_p0, const Point& trngl_p1, const Point& trngl_p2 );

Point  linesClosestPoint( const Point& line0_p0, const Point& line0_p1,
                          const Point& line1_p0, const Point& line1_p1 );

double linesDistance(    const Point& line0_p0, const Point& line0_p1,
                         const Point& line1_p0, const Point& line1_p1 );

double segmentsDistance( const Point& segm0_p0, const Point& segm0_p1,
                         const Point& segm1_p0, const Point& segm1_p1 );


// CPA - Closest Point of Approach.
double cpaTime(     const Point& start0, const Vec& vel0,
                    const Point& start1, const Vec& vel1 );

double cpaDistance( const Point& start0, const Vec& vel0,
                    const Point& start1, const Vec& vel1 );

} // namespace spatalgs
} // namespace tva
