#include "SpatialAlgorithms.h"

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


Point3 spatialalgs::project(const Point3& point, const Line3& line)
{
	return line.points[0] + (point - line.points[0]).project(line.points[1] - line.points[0]);
}

const bool spatialalgs::project(Point3& out, const Point3& point, const Segment3& segment)
{
	Vec3 res = segment.points[0] + (point - segment.points[0]).project(segment.points[1] - segment.points[0]);

	const double EPS = (segment.points[1] - segment.points[0]).magnitude() * 1e-6;
	if (INSIDE_RECTANGLE(
			segment.points[0],
			segment.points[1],
			res))
	{
		out = res;
		return true;
	}

	return false;
}


const bool spatialalgs::rayIntersectPlane(Point3& out_intersectPoint, const Ray3& ray, const Point3& pl_p0, const Point3& pl_p1, const Point3& pl_p2)
{
	Vec3 edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

	Vec3 pvec = Vec3::crossProduct(ray.direction, edges[1]);
	double det = Vec3::dotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	Vec3 tvec = ray.origin - pl_p0;
	Vec3 qvec = Vec3::crossProduct(tvec, edges[0]);

	double t = Vec3::dotProduct(edges[1], qvec) / det;
	out_intersectPoint = ray.origin + t * ray.direction;
	return t > 0.0;
}

const bool spatialalgs::isRayIntersectTriangle(const Ray3& ray, const Vec3& tr_p0, const Vec3& tr_p1, const Vec3& tr_p2)
{
	Vec3 edges[2]{ tr_p1 - tr_p0, tr_p2 - tr_p0 };

	Vec3 pvec = Vec3::crossProduct(ray.direction, edges[1]);
	double det = Vec3::dotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	double inv_det = 1.0 / det;

	Vec3 tvec = ray.origin - tr_p0;
	double u = Vec3::dotProduct(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	Vec3 qvec = Vec3::crossProduct(tvec, edges[0]);
	double v = Vec3::dotProduct(ray.direction, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	double t = Vec3::dotProduct(edges[1], qvec) * inv_det;
	return t > 0.0;
}

const bool spatialalgs::lineIntersectPlane(Vec3& out_intersectPoint, const Line3& line, const Vec3& pl_p0, const Vec3& pl_p1, const Vec3& pl_p2)
{
	const Vec3 dir = line.points[1] - line.points[0];

	Vec3 edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

	Vec3 pvec = Vec3::crossProduct(dir, edges[1]);
	double det = Vec3::dotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	Vec3 tvec = line.points[0] - pl_p0;
	Vec3 qvec = Vec3::crossProduct(tvec, edges[0]);

	double t = Vec3::dotProduct(edges[1], qvec) / det;
	out_intersectPoint = line.points[0] + t * dir;
	return true;
}

const bool spatialalgs::isSegmentIntersectTriangle(const Segment3& segment, const Vec3& tr_p0, const Vec3& tr_p1, const Vec3& tr_p2)
{
	const Vec3 dir = segment.points[1] - segment.points[0];

	Vec3 edges[2]{ tr_p1 - tr_p0, tr_p2 - tr_p0 };

	Vec3 pvec = Vec3::crossProduct(dir, edges[1]);
	double det = Vec3::dotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	double inv_det = 1.0 / det;

	Vec3 tvec = segment.points[0] - tr_p0;
	double u = Vec3::dotProduct(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	Vec3 qvec = Vec3::crossProduct(tvec, edges[0]);
	double v = Vec3::dotProduct(dir, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	double t = Vec3::dotProduct(edges[1], qvec) * inv_det;
	return t < 1.0 && t > 0.0;
}

const bool spatialalgs::segmentIntersectPlane(Vec3& out_intersectPoint, const Segment3& segment, const Vec3& pl_p0, const Vec3& pl_p1, const Vec3& pl_p2)
{
	const Vec3 dir = segment.points[1] - segment.points[0];

	Vec3 edges[2]{ pl_p1 - pl_p0, pl_p2 - pl_p0 };

	Vec3 pvec = Vec3::crossProduct(dir, edges[1]);
	double det = Vec3::dotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	Vec3 tvec = segment.points[0] - pl_p0;
	Vec3 qvec = Vec3::crossProduct(tvec, edges[0]);

	double t = Vec3::dotProduct(edges[1], qvec) / det;
	out_intersectPoint = segment.points[0] + t * dir;
	return t < 1.0 && t > 0.0;
}


const double spatialalgs::linesDistance(const Line3& line0, const Line3& line1)
{
	Vec3 u = line0.points[1] - line0.points[0];
	Vec3 v = line1.points[1] - line1.points[0];
	Vec3 w = line0.points[0] - line1.points[0];
	double a = Vec3::dotProduct(u, u); // Always >= 0
	double b = Vec3::dotProduct(u, v);
	double c = Vec3::dotProduct(v, v); // Always >= 0
	double d = Vec3::dotProduct(u, w);
	double e = Vec3::dotProduct(v, w);
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
		sc = DET(b, c, d, e) / D;
		tc = DET(a, b, d, e) / D;
	}

	// Get the difference of the two closest points
	Vec3 diff_p = w + (sc * u) - (tc * v);  // =  L1(sc) - L2(tc)

	return diff_p.magnitude();   // Return the closest distance
}

const double spatialalgs::segmentsDistance(const Segment3& segment0, const Segment3& segment1)
{
	Vec3 u = segment0.points[1] - segment0.points[0];
	Vec3 v = segment1.points[1] - segment1.points[0];
	Vec3 w = segment0.points[0] - segment1.points[0];
	double a = Vec3::dotProduct(u, u); // Always >= 0
	double b = Vec3::dotProduct(u, v);
	double c = Vec3::dotProduct(v, v); // Always >= 0
	double d = Vec3::dotProduct(u, w);
	double e = Vec3::dotProduct(v, w);
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

		if (sN < 0.0) // sc < 0 => the s=0 edge is visible
		{
			sN = 0.0;
			tN = e;
			tD = c;
		}
		else if (sN > sD) // sc > 1  => the s=1 edge is visible
		{
			sN = sD;
			tN = e + b;
			tD = c;
		}
	}

	if (tN < 0.0) // tc < 0 => the t=0 edge is visible
	{
		tN = 0.0;

		// Recompute sc for this edge
		if (-d < 0.0)
			sN = 0.0;
		else if (-d > a)
			sN = sD;
		else
		{
			sN = -d;
			sD = a;
		}
	}
	else if (tN > tD) // tc > 1  => the t=1 edge is visible
	{
		tN = tD;

		// Recompute sc for this edge
		if (b - d < 0.0)
			sN = 0;
		else if (b - d > a)
			sN = sD;
		else
		{
			sN = b - d;
			sD = a;
		}
	}

	// Finally do the division to get sc and tc
	sc = (abs(sN) < 1e-6 ? 0.0 : sN / sD);
	tc = (abs(tN) < 1e-6 ? 0.0 : tN / tD);

	// Get the difference of the two closest points
	Vec3 diff_p = w + (sc * u) - (tc * v); // =  S1(sc) - S2(tc)

	return diff_p.magnitude(); // Return the closest distance
}


const double spatialalgs::cpaTime(const Track3& track0, const Track3& track1)
{
	Vec3 dv = track0.velocity - track1.velocity;
	double dv2 = Vec3::dotProduct(dv, dv);
	if (dv2 < 1e-6) // The tracks are almost parallel.
		return 0.0;

	Vec3 w0 = track0.origin - track1.origin;
	return -Vec3::dotProduct(w0, dv) / dv2;
}

const double spatialalgs::cpaDistance(const Track3& track0, const Track3& track1)
{
	double cpa_time = cpaTime({ track0.origin, track0.velocity }, { track1.origin, track1.velocity });
	Vec3 p0 = track0.origin + (cpa_time * track0.velocity);
	Vec3 p1 = track1.origin + (cpa_time * track1.velocity);
	return (p1 - p0).magnitude();
}
