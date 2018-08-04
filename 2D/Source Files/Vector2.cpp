#include "Vector2.h"


const double Vector2::GetCoordinate(const int& axisIndex)
{
	return _coors[axisIndex];
}

Vector2& Vector2::operator=(const Vector2& vec)
{
	_coors[0] = vec._coors[0];
	_coors[1] = vec._coors[1];

	return *this;
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