#pragma once
#include <cmath>
#include <algorithm>
#include "Vec3.h"
#include "Segment3.h"
#include "Ray3.h"
#include "Track3.h"


namespace spatialalgs
{
	      Point3 project(const Point3& point, const Line3& line);
	const bool   project(Point3& out, const Point3& point, const Segment3& segment);

	const bool isRayIntersectPlane(const Ray3& ray, const Vec3& pl_p0, const Vec3& pl_p1, const Vec3& pl_p2);
	const bool rayIntersectPlane(Point3& out_intersectPoint, const Ray3& ray, const Point3& pl_p0, const Point3& pl_p1, const Point3& pl_p2);
	const bool isRayIntersectTriangle(const Ray3& ray, const Vec3& tr_p0, const Vec3& tr_p1, const Vec3& tr_p2);
	const bool rayIntersectTriangle(Point3& out_intersectPoint, const Ray3& ray, const Point3& tr_p0, const Point3& tr_p1, const Point3& tr_p2);
	const bool isLineIntersectPlane(const Ray3& ray, const Vec3& pl_p0, const Vec3& pl_p1, const Vec3& pl_p2);
	const bool lineIntersectPlane(Vec3& out_intersectPoint, const Line3& line, const Vec3& pl_p0, const Vec3& pl_p1, const Vec3& pl_p2);
	const bool isSegmentIntersectTriangle(const Segment3& segment, const Vec3& tr_p0, const Vec3& tr_p1, const Vec3& tr_p2);
	const bool segmentIntersectTriangle(Vec3& out_intersectPoint, const Segment3& segment, const Vec3& tr_p0, const Vec3& tr_p1, const Vec3& tr_p2);
	const bool isSegmentIntersectPlane(const Segment3& segment, const Vec3& pl_p0, const Vec3& pl_p1, const Vec3& pl_p2);
	const bool segmentIntersectPlane(Vec3& out_intersectPoint, const Segment3& segment, const Vec3& pl_p0, const Vec3& pl_p1, const Vec3& pl_p2);

	const double linesDistance(const Line3& line0, const Line3& line1);
	const double segmentsDistance(const Segment3& segment0, const Segment3& segment1);

	// CPA - Closest Point of Approach.
	const double cpaTime(const Track3& track0, const Track3& track1);
	const double cpaDistance(const Track3& track0, const Track3& track1);
}