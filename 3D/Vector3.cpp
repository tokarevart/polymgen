#include "Vector3.h"
#include <algorithm>

#define DET(a, b, c, d) \
		((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0._coors[0], corner1._coors[0], point._coors[0]) && \
		 BETWEEN(corner0._coors[1], corner1._coors[1], point._coors[1]) && \
		 BETWEEN(corner0._coors[2], corner1._coors[2], point._coors[2]))

#define PI 3.141592653589793

const double PI_DIV_180 = 3.141592653589793 / 180.0;

double Sign(double val)
{
	return val > 0.0 ? 1.0 : -1.0;
}

const double Vector3::DotProduct(
	const Vector3& vec0, 
	const Vector3& vec1)
{
	return vec0._coors[0] * vec1._coors[0] + vec0._coors[1] * vec1._coors[1] + vec0._coors[2] * vec1._coors[2];
}

Vector3 Vector3::CrossProduct(
	const Vector3& vec0, 
	const Vector3& vec1)
{
	return Vector3(
		vec0._coors[1] * vec1._coors[2] - vec0._coors[2] * vec1._coors[1],
		vec0._coors[2] * vec1._coors[0] - vec0._coors[0] * vec1._coors[2],
		vec0._coors[0] * vec1._coors[1] - vec0._coors[1] * vec1._coors[0]);
}

const double Vector3::MixedProduct(const Vector3 &vec0, const Vector3 &vec1, const Vector3 &vec2)
{
	return DotProduct(CrossProduct(vec0, vec1), vec2);
}

const double Vector3::Cos(
	const Vector3& vec0, 
	const Vector3& vec1)
{
	return Vector3::DotProduct(vec0, vec1) / sqrt(vec0.SqrMagnitude() * vec1.SqrMagnitude());
}

Vector3 Vector3::Project(
	const Vector3& point,
	const Vector3& line_p0, 
	const Vector3& line_p1)
{
	return line_p0 + Vector3(point - line_p0).Project(line_p1 - line_p0);
}

const bool Vector3::Project(
	Vector3& out, 
	const Vector3& point, 
	const Vector3& segm_p0, 
	const Vector3& segm_p1)
{
	Vector3 res = segm_p0 + (point - segm_p0).Project(segm_p1 - segm_p0);

	const double EPS = (segm_p1 - segm_p0).Magnitude() * 1e-6;
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

const bool Vector3::RayIntersectTriangle(
	const Vector3 &origin, 
	const Vector3 &dir, 
	const Vector3 &tr_p0, 
	const Vector3 &tr_p1, 
	const Vector3 &tr_p2)
{
	Vector3 edges[2] 
		{ tr_p1 - tr_p0, 
		  tr_p2 - tr_p0 };

	Vector3 pvec = CrossProduct(dir, edges[1]);
	double det = DotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	double inv_det = 1.0 / det;

	Vector3 tvec = origin - tr_p0;
	double u = DotProduct(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	Vector3 qvec = CrossProduct(tvec, edges[0]);
	double v = DotProduct(dir, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	double t = DotProduct(edges[1], qvec) * inv_det;
	return t > 0.0 ? true : false;
}

const bool Vector3::SegmentIntersectTriangle(
	const Vector3 &segm_p0, 
	const Vector3 &segm_p1,
	const Vector3 &tr_p0, 
	const Vector3 &tr_p1, 
	const Vector3 &tr_p2)
{
	const Vector3 dir = segm_p1 - segm_p0;

	Vector3 edges[2]
	{ tr_p1 - tr_p0,
		tr_p2 - tr_p0 };

	Vector3 pvec = CrossProduct(dir, edges[1]);
	double det = DotProduct(edges[0], pvec);
	if (det < 1e-6 && det > -1e-6)
		return false;

	double inv_det = 1.0 / det;

	Vector3 tvec = segm_p0 - tr_p0;
	double u = DotProduct(tvec, pvec) * inv_det;
	if (u < 0.0 || u > 1.0)
		return false;

	Vector3 qvec = CrossProduct(tvec, edges[0]);
	double v = DotProduct(dir, qvec) * inv_det;
	if (v < 0.0 || u + v > 1.0)
		return false;

	double t = DotProduct(edges[1], qvec) * inv_det;
	return t < 1.0 && t > 0.0;
}

// I didn't tested it.
const double Vector3::LinesDistance(const Vector3 &l0_p0, const Vector3 &l0_p1, const Vector3 &l1_p0, const Vector3 &l1_p1)
{
	Vector3 u = l0_p1 - l0_p0;
	Vector3 v = l1_p1 - l1_p0;
	Vector3 w = l0_p0 - l1_p0;
	double a = DotProduct(u, u); // Always >= 0
	double b = DotProduct(u, v);
	double c = DotProduct(v, v); // Always >= 0
	double d = DotProduct(u, w);
	double e = DotProduct(v, w);
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
	Vector3 diff_p = w + (sc * u) - (tc * v);  // =  L1(sc) - L2(tc)

	return diff_p.Magnitude();   // Return the closest distance
}

// I didn't tested it.
const double Vector3::SegmentsDistance(const Vector3 &s0_p0, const Vector3 &s0_p1, const Vector3 &s1_p0, const Vector3 &s1_p1)
{
	Vector3 u = s0_p1 - s0_p0;
	Vector3 v = s1_p1 - s1_p0;
	Vector3 w = s0_p0 - s1_p0;
	double a = DotProduct(u, u); // Always >= 0
	double b = DotProduct(u, v);
	double c = DotProduct(v, v); // Always >= 0
	double d = DotProduct(u, w);
	double e = DotProduct(v, w);
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
	Vector3 diff_p = w + (sc * u) - (tc * v); // =  S1(sc) - S2(tc)

	return diff_p.Magnitude(); // Return the closest distance
}

const double Vector3::Magnitude() const
{
	return sqrt(SqrMagnitude());
}

const double Vector3::SqrMagnitude() const
{
	return DotProduct(*this, *this);
}

Vector3& Vector3::Normalize()
{
	double inv_magn = 1.0 / Magnitude();
	_coors[0] *= inv_magn;
	_coors[1] *= inv_magn;
	_coors[2] *= inv_magn;
	return *this;
}

Vector3& Vector3::Project(const Vector3& vec)
{
	*this = (DotProduct(*this, vec) / vec.SqrMagnitude()) * vec;
	return *this;
}

Vector3& Vector3::Project(
	const Vector3& plane_v0, 
	const Vector3& plane_v1)
{
	*this -= Vector3(*this).Project(CrossProduct(plane_v0, plane_v1));
	return *this;
}

double Vector3::DistanceToLine(const Vector3 &line_p0, const Vector3 &line_p1) const
{
	return (Project(*this, line_p0, line_p1) - *this).Magnitude();
}

double Vector3::DistanceToSegment(const Vector3 &segm_p0, const Vector3 &segm_p1) const
{
	Vector3 proj = Project(*this, segm_p0, segm_p1);
	double sqr_magns[2];

	const double EPS = (segm_p1 - segm_p0).Magnitude() * 1e-6;
	if (INSIDE_RECTANGLE(segm_p0, segm_p1, proj))
		return (proj - *this).Magnitude();

	else if ((sqr_magns[0] = (segm_p0 - *this).SqrMagnitude()) < (sqr_magns[1] = (segm_p1 - *this).SqrMagnitude()))
		return sqrt(sqr_magns[0]);

	else
		return sqrt(sqr_magns[1]);
}

Vector3& Vector3::operator=(const Vector3& vec)
{
	if (this == &vec)
		return *this;

	_coors[0] = vec._coors[0];
	_coors[1] = vec._coors[1];
	_coors[2] = vec._coors[2];
	return *this;
}

double& Vector3::operator[](const int& axisIndex)
{
	return _coors[axisIndex];
}

Vector3 Vector3::operator+() const
{
	return Vector3(*this);
}

Vector3 Vector3::operator-() const
{
	return Vector3(-_coors[0], -_coors[1], -_coors[2]);
}

Vector3 Vector3::operator+(const Vector3& right) const
{
	return Vector3(_coors[0] + right._coors[0], _coors[1] + right._coors[1], _coors[2] + right._coors[2]);
}

Vector3 Vector3::operator-(const Vector3& right) const
{
	return Vector3(_coors[0] - right._coors[0], _coors[1] - right._coors[1], _coors[2] - right._coors[2]);
}

Vector3& Vector3::operator+=(const Vector3& right)
{
	_coors[0] += right._coors[0];
	_coors[1] += right._coors[1];
	_coors[2] += right._coors[2];
	return *this;
}

Vector3& Vector3::operator-=(const Vector3& right)
{
	_coors[0] -= right._coors[0];
	_coors[1] -= right._coors[1];
	_coors[2] -= right._coors[2];
	return *this;
}

Vector3& Vector3::operator*=(const double& scalar)
{
	_coors[0] *= scalar;
	_coors[1] *= scalar;
	_coors[2] *= scalar;
	return *this;
}

Vector3& Vector3::operator/=(const double& scalar)
{
	_coors[0] /= scalar;
	_coors[1] /= scalar;
	_coors[2] /= scalar;
	return *this;
}

const Vector3 operator*(
	const Vector3& vec, 
	const double& scalar)
{
	return Vector3(
		vec._coors[0] * scalar, 
		vec._coors[1] * scalar, 
		vec._coors[2] * scalar);
}

const Vector3 operator*(
	const double& scalar, 
	const Vector3& vec)
{
	return Vector3(
		scalar * vec._coors[0], 
		scalar * vec._coors[1], 
		scalar * vec._coors[2]);
}

const Vector3 operator/(
	const Vector3& vec, 
	const double& scalar)
{
	double inv_scalar = 1.0 / scalar;
	return Vector3(
		vec._coors[0] * inv_scalar, 
		vec._coors[1] * inv_scalar, 
		vec._coors[2] * inv_scalar);
}

Vector3::Vector3() 
{
	_coors[0] = 0.0;
	_coors[1] = 0.0;
	_coors[2] = 0.0;
}

Vector3::Vector3(const Vector3& vec)
{
	*this = vec;
}

Vector3::Vector3(
	const double& coor0, 
	const double& coor1, 
	const double& coor2)
{
	_coors[0] = coor0;
	_coors[1] = coor1;
	_coors[2] = coor2;
}

Vector3::~Vector3() {}