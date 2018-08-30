#include "Vector3.h"
#include <algorithm>

#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0._coors[0], corner1._coors[0], point._coors[0]) && \
		 BETWEEN(corner0._coors[1], corner1._coors[1], point._coors[1]) && \
		 BETWEEN(corner0._coors[2], corner1._coors[2], point._coors[2]))

const double PI_DIV_180 = 3.141592653589793 / 180.0;


double Sign(double val)
{
	return val > 0.0 ? 1.0 : -1.0;
}

const double Vector3::DotProduct(const Vector3& vec0, const Vector3& vec1)
{
	return vec0._coors[0] * vec1._coors[0] + vec0._coors[1] * vec1._coors[1] + vec0._coors[2] * vec1._coors[2];
}

Vector3 Vector3::CrossProduct(const Vector3& vec0, const Vector3& vec1)
{
	return Vector3(
		vec0._coors[1] * vec1._coors[2] - vec0._coors[2] * vec1._coors[1],
		vec0._coors[2] * vec1._coors[0] - vec0._coors[0] * vec1._coors[2],
		vec0._coors[0] * vec1._coors[1] - vec0._coors[1] * vec1._coors[0]);
}

const double Vector3::Cos(const Vector3& vec0, const Vector3& vec1)
{
	return Vector3::DotProduct(vec0, vec1) / sqrt(vec0.SqrMagnitude() * vec1.SqrMagnitude());
}

Vector3 Vector3::Project(const Vector3& point, const Vector3& linePoint0, const Vector3& linePoint1)
{
	return linePoint0 + Vector3(point - linePoint0).Project(linePoint1 - linePoint0);
}

const bool Vector3::Project(Vector3& out, const Vector3& point, const Vector3& segmPoint0, const Vector3& segmPoint1)
{
	Vector3 res(point);
	res.Project(segmPoint0, segmPoint1);

	const double EPS = (segmPoint1 - segmPoint0).Magnitude() * 1e-6;
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
	double buf = 1.0 / Magnitude();
	_coors[0] *= buf;
	_coors[1] *= buf;
	_coors[2] *= buf;

	return *this;
}

Vector3& Vector3::Project(const Vector3& vec)
{
	Vector3 res = (DotProduct(*this, vec) / vec.SqrMagnitude()) * vec;
	_coors[0] = res._coors[0];
	_coors[1] = res._coors[1];

	return *this;
}

Vector3& Vector3::Project(const Vector3& planeVec0, const Vector3& planeVec1)
{
	Vector3 cross_buf = CrossProduct(planeVec0, planeVec1);
	*this = *this - Vector3(*this).Project(cross_buf);

	return *this;
}

Vector3& Vector3::operator=(const Vector3& vec)
{
	if (this == &vec)
	{
		return *this;
	}

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

const Vector3 operator*(const Vector3& vec, const double& scalar)
{
	return Vector3(vec._coors[0] * scalar, vec._coors[1] * scalar, vec._coors[2] * scalar);
}

const Vector3 operator*(const double& scalar, const Vector3& vec)
{
	return Vector3(scalar * vec._coors[0], scalar * vec._coors[1], scalar * vec._coors[2]);
}

const Vector3 operator/(const Vector3& vec, const double& scalar)
{
	double inv_scalar = 1.0 / scalar;
	return Vector3(vec._coors[0] * inv_scalar, vec._coors[1] * inv_scalar, vec._coors[2] * inv_scalar);
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

Vector3::Vector3(const double& coor0, const double& coor1, const double& coor2)
{
	_coors[0] = coor0;
	_coors[1] = coor1;
	_coors[2] = coor2;
}

Vector3::~Vector3() {}