#include "Vector2.h"
#include <algorithm>


#define DET(a, b, c, d) \
		(a * d - b * c)

#define EPS 1e-10 // Must depend on shell edges lengths
#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0._coors[0], corner1._coors[0], point._coors[0]) && \
		 BETWEEN(corner0._coors[1], corner1._coors[1], point._coors[1]))

#define NOT_NEAR_SEGMS(segm0p0, segm0p1, segm1p0, segm1p1) \
		(segm0p0._coors[0] < segm1p0._coors[0]  && \
		 segm0p1._coors[0] < segm1p0._coors[0] && \
		 segm0p0._coors[0] < segm1p1._coors[0] && \
		 segm0p1._coors[0] < segm1p1._coors[0]) || \
		\
		(segm0p0._coors[1] < segm1p0._coors[1] && \
		 segm0p1._coors[1] < segm1p0._coors[1] && \
		 segm0p0._coors[1] < segm1p1._coors[1] && \
		 segm0p1._coors[1] < segm1p1._coors[1]) || \
		\
		(segm0p0._coors[0] > segm1p0._coors[0]  && \
		 segm0p1._coors[0] > segm1p0._coors[0] && \
		 segm0p0._coors[0] > segm1p1._coors[0] && \
		 segm0p1._coors[0] > segm1p1._coors[0]) || \
		\
		(segm0p0._coors[1] > segm1p0._coors[1] && \
		 segm0p1._coors[1] > segm1p0._coors[1] && \
		 segm0p0._coors[1] > segm1p1._coors[1] && \
		 segm0p1._coors[1] > segm1p1._coors[1])

const double PI_DIV_180 = 3.141592653589793 / 180.0;


const double Vector2::DotProduct(const Vector2& vec0, const Vector2& vec1)
{
	return vec0._coors[0] * vec1._coors[0] + vec0._coors[1] * vec1._coors[1];
}

const double Vector2::CrossProductMagnitude(const Vector2& vec0, const Vector2& vec1)
{
	return vec0._coors[0] * vec1._coors[1] - vec0._coors[1] * vec1._coors[0];
}

const double Vector2::Cos(const Vector2& vec0, const Vector2& vec1)
{
	return Vector2::DotProduct(vec0, vec1) / sqrt(vec0.SqrMagnitude() * vec1.SqrMagnitude());
}

Vector2 Vector2::Project(const Vector2& point, const Vector2& linePoint0, const Vector2& linePoint1)
{
	return Vector2(point).Project(linePoint0, linePoint1);
}

const bool Vector2::Project(Vector2& out, const Vector2& point, const Vector2& segmPoint0, const Vector2& segmPoint1)
{
	Vector2 res(point);
	res.Project(segmPoint0, segmPoint1);

	if (INSIDE_RECTANGLE(
			segmPoint0,
			segmPoint1,
			res))
	{
		out = res;
		return true;
	}

	return false;
}

Vector2 Vector2::LinesIntersection(const Vector2& line0p0, const Vector2& line0p1, const Vector2& line1p0, const Vector2& line1p1)
{
	double A1 = line0p0._coors[1] - line0p1._coors[1];
	double B1 = line0p1._coors[0] - line0p0._coors[0];
	double C1 = -A1 * line0p0._coors[0] - B1 * line0p0._coors[1];
	double A2 = line1p0._coors[1] - line1p1._coors[1];
	double B2 = line1p1._coors[0] - line1p0._coors[0];
	double C2 = -A2 * line1p0._coors[0] - B2 * line1p0._coors[1];

	double minus_inv_zn = -1.0 / DET(A1, B1, A2, B2);

	double x = DET(C1, B1, C2, B2) * minus_inv_zn;
	double y = DET(A1, C1, A2, C2) * minus_inv_zn;
	return Vector2(x, y);
}

const bool Vector2::SegmentsIntersection(const Vector2& segm0p0, const Vector2& segm0p1, const Vector2& segm1p0, const Vector2& segm1p1)
{
	if (NOT_NEAR_SEGMS(
			segm0p0,
			segm0p1,
			segm1p0,
			segm1p1))
	{
		return false;
	}

	Vector2 intersect_point =
		LinesIntersection(
			segm0p0,
			segm0p1,
			segm1p0,
			segm1p1);

	if (INSIDE_RECTANGLE(
			segm0p0,
			segm0p1,
			intersect_point) &&
		INSIDE_RECTANGLE(
			segm1p0,
			segm1p1,
			intersect_point))
	{
		return true;
	}

	return false;
}

const bool Vector2::SegmentsIntersection(Vector2& out, const Vector2& segm0p0, const Vector2& segm0p1, const Vector2& segm1p0, const Vector2& segm1p1)
{
	if (NOT_NEAR_SEGMS(
			segm0p0,
			segm0p1,
			segm1p0,
			segm1p1))
	{
		return false;
	}

	Vector2 intersect_point =
		LinesIntersection(
			segm0p0,
			segm0p1,
			segm1p0,
			segm1p1);

	if (INSIDE_RECTANGLE(
			segm0p0,
			segm0p1,
			intersect_point) &&
		INSIDE_RECTANGLE(
			segm1p0,
			segm1p1,
			intersect_point))
	{
		out = intersect_point;
		return true;
	}

	return false;
}

const double Vector2::Magnitude() const
{
	return sqrt(SqrMagnitude());
}

const double Vector2::SqrMagnitude() const
{
	return DotProduct(*this, *this);
}

Vector2& Vector2::Normalize()
{
	double buf = 1.0 / Magnitude();
	_coors[0] *= buf;
	_coors[1] *= buf;

	return *this;
}

Vector2& Vector2::Rotate(const double& angle, const AngleUnit& unit)
{
	double buf = _coors[0];
	switch (unit)
	{
	case Radian:
		_coors[0] = cos(angle) * _coors[0] - sin(angle) * _coors[1];
		_coors[1] = sin(angle) * buf       + cos(angle) * _coors[1];
		break;

	case Degree:
		double angle_rad = angle * PI_DIV_180;
		_coors[0] = cos(angle_rad) * _coors[0] - sin(angle_rad) * _coors[1];
		_coors[1] = sin(angle_rad) * buf       + cos(angle_rad) * _coors[1];
		break;

	default:
		break;
	}

	return *this;
}

Vector2& Vector2::Mirror(const Vector2& vec)
{
	double angle = acos(Cos(*this, vec));

	return CrossProductMagnitude(vec, *this) > 0 ? 
		Rotate(-2.0 * angle, Radian) :
		Rotate(2.0 * angle, Radian);
}

Vector2& Vector2::Project(const Vector2& vec)
{
	Vector2 res = (DotProduct(*this, vec) / vec.SqrMagnitude()) * vec;
	_coors[0] = res._coors[0];
	_coors[1] = res._coors[1];

	return *this;
}

Vector2& Vector2::Project(const Vector2& linePoint0, const Vector2& linePoint1)
{
	*this = linePoint0 + Vector2(*this - linePoint0).Project(linePoint1 - linePoint0);

	return *this;
}

Vector2& Vector2::operator=(const Vector2& vec)
{
	if (this == &vec)
	{
		return *this;
	}

	_coors[0] = vec._coors[0];
	_coors[1] = vec._coors[1];

	return *this;
}

double& Vector2::operator[](const int& axisIndex)
{
	return _coors[axisIndex];
}

Vector2 Vector2::operator+() const
{
	return Vector2(*this);
}

Vector2 Vector2::operator-() const
{
	return Vector2(-_coors[0], -_coors[1]);
}

Vector2 Vector2::operator+(const Vector2& right) const
{
	return Vector2(_coors[0] + right._coors[0], _coors[1] + right._coors[1]);
}

Vector2 Vector2::operator-(const Vector2& right) const
{
	return Vector2(_coors[0] - right._coors[0], _coors[1] - right._coors[1]);
}

Vector2& Vector2::operator+=(const Vector2& right)
{
	_coors[0] += right._coors[0];
	_coors[1] += right._coors[1];

	return *this;
}

Vector2& Vector2::operator-=(const Vector2& right)
{
	_coors[0] -= right._coors[0];
	_coors[1] -= right._coors[1];

	return *this;
}

Vector2& Vector2::operator*=(const double& scalar)
{
	_coors[0] *= scalar;
	_coors[1] *= scalar;

	return *this;
}

const Vector2 operator*(const Vector2& vec, const double& scalar)
{
	return Vector2(vec._coors[0] * scalar, vec._coors[1] * scalar);
}

const Vector2 operator*(const double& scalar, const Vector2& vec)
{
	return Vector2(scalar * vec._coors[0], scalar * vec._coors[1]);
}

const Vector2 operator/(const Vector2 & vec, const double & scalar)
{
	double inv_scalar = 1.0 / scalar;
	return Vector2(vec._coors[0] * inv_scalar, vec._coors[1] * inv_scalar);
}

Vector2::Vector2() 
{
	_coors[0] = 0.0;
	_coors[1] = 0.0;
}

Vector2::Vector2(const Vector2& vec)
{
	*this = vec;
}

Vector2::Vector2(const double& coor0, const double& coor1)
{
	_coors[0] = coor0;
	_coors[1] = coor1;
}

Vector2::~Vector2() {}