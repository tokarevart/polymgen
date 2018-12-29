#include "Vec3.h"
#include <algorithm>

#define DET(a, b, c, d) \
		((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0._coors[0], corner1._coors[0], point._coors[0]) && \
		 BETWEEN(corner0._coors[1], corner1._coors[1], point._coors[1]) && \
		 BETWEEN(corner0._coors[2], corner1._coors[2], point._coors[2]))

#define PI         3.141592653589793
#define PI_DIV_180 0.017453292519943295


const double Vec3::dotProduct(const Vec3& vec0, const Vec3& vec1)
{
	return vec0._coors[0] * vec1._coors[0] + vec0._coors[1] * vec1._coors[1] + vec0._coors[2] * vec1._coors[2];
}

Vec3 Vec3::crossProduct(const Vec3& vec0, const Vec3& vec1)
{
	return Vec3(
		vec0._coors[1] * vec1._coors[2] - vec0._coors[2] * vec1._coors[1],
		vec0._coors[2] * vec1._coors[0] - vec0._coors[0] * vec1._coors[2],
		vec0._coors[0] * vec1._coors[1] - vec0._coors[1] * vec1._coors[0]);
}

const double Vec3::mixedProduct(const Vec3 &vec0, const Vec3 &vec1, const Vec3 &vec2)
{
	return dotProduct(crossProduct(vec0, vec1), vec2);
}

const double Vec3::cos(const Vec3& vec0, const Vec3& vec1)
{
	return Vec3::dotProduct(vec0, vec1) / sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
}

Vec3 Vec3::project(const Vec3& point, const Vec3& line_p0, const Vec3& line_p1)
{
	return line_p0 + Vec3(point - line_p0).project(line_p1 - line_p0);
}

const bool Vec3::project(
	Vec3& out, 
	const Vec3& point, 
	const Vec3& segm_p0, 
	const Vec3& segm_p1)
{
	Vec3 res = segm_p0 + (point - segm_p0).project(segm_p1 - segm_p0);

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

const bool Vec3::rayIntersectTriangle(
	const Vec3 &origin, 
	const Vec3 &dir, 
	const Vec3 &tr_p0, 
	const Vec3 &tr_p1, 
	const Vec3 &tr_p2)
{
	Vec3 edges[2] 
		{ tr_p1 - tr_p0, 
		  tr_p2 - tr_p0 };

	Vec3 pvec = crossProduct(dir, edges[1]);
	double det = dotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	double inv_det = 1.0 / det;

	Vec3 tvec = origin - tr_p0;
	double u = dotProduct(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	Vec3 qvec = crossProduct(tvec, edges[0]);
	double v = dotProduct(dir, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	double t = dotProduct(edges[1], qvec) * inv_det;
	return t > 0.0 ? true : false;
}

const bool Vec3::segmentIntersectTriangle(
	const Vec3 &segm_p0, 
	const Vec3 &segm_p1,
	const Vec3 &tr_p0, 
	const Vec3 &tr_p1, 
	const Vec3 &tr_p2)
{
	const Vec3 dir = segm_p1 - segm_p0;

	Vec3 edges[2]
	{ tr_p1 - tr_p0,
		tr_p2 - tr_p0 };

	Vec3 pvec = crossProduct(dir, edges[1]);
	double det = dotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	double inv_det = 1.0 / det;

	Vec3 tvec = segm_p0 - tr_p0;
	double u = dotProduct(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	Vec3 qvec = crossProduct(tvec, edges[0]);
	double v = dotProduct(dir, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	double t = dotProduct(edges[1], qvec) * inv_det;
	return t < 1.0 && t > 0.0;
}

// I didn't tested it.
const double Vec3::linesDistance(const Vec3 &l0_p0, const Vec3 &l0_p1, const Vec3 &l1_p0, const Vec3 &l1_p1)
{
	Vec3 u = l0_p1 - l0_p0;
	Vec3 v = l1_p1 - l1_p0;
	Vec3 w = l0_p0 - l1_p0;
	double a = dotProduct(u, u); // Always >= 0
	double b = dotProduct(u, v);
	double c = dotProduct(v, v); // Always >= 0
	double d = dotProduct(u, w);
	double e = dotProduct(v, w);
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

// I didn't tested it.
const double Vec3::segmentsDistance(const Vec3 &s0_p0, const Vec3 &s0_p1, const Vec3 &s1_p0, const Vec3 &s1_p1)
{
	Vec3 u = s0_p1 - s0_p0;
	Vec3 v = s1_p1 - s1_p0;
	Vec3 w = s0_p0 - s1_p0;
	double a = dotProduct(u, u); // Always >= 0
	double b = dotProduct(u, v);
	double c = dotProduct(v, v); // Always >= 0
	double d = dotProduct(u, w);
	double e = dotProduct(v, w);
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

const double Vec3::cpaTime(const Vec3& start0, const Vec3& velocity0, const Vec3& start1, const Vec3& velocity1)
{
	Vec3 dv = velocity0 - velocity1;
	double dv2 = dotProduct(dv, dv);
	if (dv2 < 1e-6) // The tracks are almost parallel.
		return 0.0;

	Vec3 w0 = start0 - start1;
	return -dotProduct(w0, dv) / dv2;
}

const double Vec3::cpaDistance(const Vec3& start0, const Vec3& velocity0, const Vec3& start1, const Vec3& velocity1)
{
	double cpa_time = cpaTime(start0, velocity0, start1, velocity1);
	Vec3 p0 = start0 + (cpa_time * velocity0);
	Vec3 p1 = start1 + (cpa_time * velocity1);
	return (p1 - p0).magnitude();
}

const double Vec3::magnitude() const
{
	return sqrt(sqrMagnitude());
}

const double Vec3::sqrMagnitude() const
{
	return dotProduct(*this, *this);
}

Vec3& Vec3::normalize()
{
	double inv_magn = 1.0 / magnitude();
	_coors[0] *= inv_magn;
	_coors[1] *= inv_magn;
	_coors[2] *= inv_magn;
	return *this;
}

Vec3& Vec3::project(const Vec3& vec)
{
	*this = (dotProduct(*this, vec) / vec.sqrMagnitude()) * vec;
	return *this;
}

Vec3& Vec3::project(const Vec3& plane_v0, const Vec3& plane_v1)
{
	*this -= Vec3(*this).project(crossProduct(plane_v0, plane_v1));
	return *this;
}

double Vec3::distanceToLine(const Vec3 &line_p0, const Vec3 &line_p1) const
{
	return (project(*this, line_p0, line_p1) - *this).magnitude();
}

double Vec3::distanceToSegment(const Vec3 &segm_p0, const Vec3 &segm_p1) const
{
	Vec3 proj = project(*this, segm_p0, segm_p1);
	double sqr_magns[2];

	const double EPS = (segm_p1 - segm_p0).magnitude() * 1e-6;
	if (INSIDE_RECTANGLE(segm_p0, segm_p1, proj))
		return (proj - *this).magnitude();

	else if ((sqr_magns[0] = (segm_p0 - *this).sqrMagnitude()) < (sqr_magns[1] = (segm_p1 - *this).sqrMagnitude()))
		return sqrt(sqr_magns[0]);

	else
		return sqrt(sqr_magns[1]);
}

Vec3& Vec3::operator=(const Vec3& vec)
{
	if (this == &vec)
		return *this;

	_coors[0] = vec._coors[0];
	_coors[1] = vec._coors[1];
	_coors[2] = vec._coors[2];
	return *this;
}

double& Vec3::operator[](const int axisIndex)
{
	return _coors[axisIndex];
}

Vec3 Vec3::operator+() const
{
	return Vec3(*this);
}

Vec3 Vec3::operator-() const
{
	return Vec3(-_coors[0], -_coors[1], -_coors[2]);
}

Vec3 Vec3::operator+(const Vec3& right) const
{
	return Vec3(_coors[0] + right._coors[0], _coors[1] + right._coors[1], _coors[2] + right._coors[2]);
}

Vec3 Vec3::operator-(const Vec3& right) const
{
	return Vec3(_coors[0] - right._coors[0], _coors[1] - right._coors[1], _coors[2] - right._coors[2]);
}

Vec3& Vec3::operator+=(const Vec3& right)
{
	_coors[0] += right._coors[0];
	_coors[1] += right._coors[1];
	_coors[2] += right._coors[2];
	return *this;
}

Vec3& Vec3::operator-=(const Vec3& right)
{
	_coors[0] -= right._coors[0];
	_coors[1] -= right._coors[1];
	_coors[2] -= right._coors[2];
	return *this;
}

Vec3& Vec3::operator*=(const double& scalar)
{
	_coors[0] *= scalar;
	_coors[1] *= scalar;
	_coors[2] *= scalar;
	return *this;
}

Vec3& Vec3::operator/=(const double& scalar)
{
	_coors[0] /= scalar;
	_coors[1] /= scalar;
	_coors[2] /= scalar;
	return *this;
}

const Vec3 operator*(const Vec3& vec, const double scalar)
{
	return Vec3(
		vec._coors[0] * scalar, 
		vec._coors[1] * scalar, 
		vec._coors[2] * scalar);
}

const Vec3 operator*(const double scalar, const Vec3& vec)
{
	return Vec3(
		scalar * vec._coors[0], 
		scalar * vec._coors[1], 
		scalar * vec._coors[2]);
}

const Vec3 operator/(const Vec3& vec, const double scalar)
{
	double inv_scalar = 1.0 / scalar;
	return Vec3(
		vec._coors[0] * inv_scalar, 
		vec._coors[1] * inv_scalar, 
		vec._coors[2] * inv_scalar);
}

Vec3::Vec3() 
{
	_coors[0] = 0.0;
	_coors[1] = 0.0;
	_coors[2] = 0.0;
}

Vec3::Vec3(const Vec3& vec)
{
	*this = vec;
}

Vec3::Vec3(const double coor0, const double coor1, const double coor2)
{
	_coors[0] = coor0;
	_coors[1] = coor1;
	_coors[2] = coor2;
}

Vec3::~Vec3() {}
