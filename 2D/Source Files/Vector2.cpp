#include "Vector2.h"


const double PI_DIV_180 = 3.141592653589793 / 180.0;

const double Vector2::GetLength()
{
	return sqrt(GetSqrLength());
}

const double Vector2::GetSqrLength()
{
	return _coors[0] * _coors[0] + _coors[1] * _coors[1];
}

Vector2& Vector2::Normalize()
{
	double buf = 1.0 / GetLength();
	_coors[0] *= buf;
	_coors[1] *= buf;

	return *this;
}

Vector2& Vector2::Rotate(const double& angle, AngleUnit unit)
{
	switch (unit)
	{
	case Radian:
		double buf = _coors[0];
		_coors[0] = cos(angle) * _coors[0] - sin(angle) * _coors[1];
		_coors[1] = sin(angle) * buf       + cos(angle) * _coors[1];
		break;

	case Degree:
		double buf = _coors[0];
		double angle_rad = angle * PI_DIV_180;
		_coors[0] = cos(angle_rad) * _coors[0] - sin(angle_rad) * _coors[1];
		_coors[1] = sin(angle_rad) * buf       + cos(angle_rad) * _coors[1];
		break;

	default:
		break;
	}

	return *this;
}

Vector2& Vector2::operator=(const Vector2& vec)
{
	_coors[0] = vec._coors[0];
	_coors[1] = vec._coors[1];

	return *this;
}

double& Vector2::operator[](const int& axisIndex)
{
	return _coors[axisIndex];
}

Vector2 Vector2::operator+()
{
	return Vector2(*this);
}

Vector2 Vector2::operator-()
{
	return Vector2(-_coors[0], -_coors[1]);
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

Vector2::~Vector2()
{
}