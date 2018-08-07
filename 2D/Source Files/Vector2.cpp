#include "Vector2.h"


const double PI_DIV_180 = 3.141592653589793 / 180.0;

const double Vector2::DotProduct(const Vector2& vec0, const Vector2& vec1)
{
	return vec0._coors[0] * vec1._coors[0] + vec0._coors[1] * vec1._coors[1];
}

const double Vector2::CrossProductMagnitude(const Vector2& vec0, const Vector2& vec1)
{
	return vec0._coors[0] * vec1._coors[1] - vec0._coors[1] * vec1._coors[0];
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
	double angle = acos(DotProduct(*this, vec) / sqrt(this->SqrMagnitude() * vec.SqrMagnitude()));

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

Vector2::Vector2() {}

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
