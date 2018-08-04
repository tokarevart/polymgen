#pragma once
#include <cmath>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

enum AngleUnit
{
	Radian,
	Degree
};

class Vector2
{
private:
	double _coors[2];

public:
	const double GetLength();
	const double GetSqrLength();
	Vector2& Normalize();
	Vector2& Rotate(const double& angle, AngleUnit unit);

	Vector2& operator=(const Vector2& vec);
	double& operator[](const int& axisIndex);
	Vector2 operator+();
	Vector2 operator-();
	Vector2 operator+(const Vector2& right);
	Vector2 operator-(const Vector2& right);
	Vector2& operator+=(const Vector2& right);
	Vector2& operator-=(const Vector2& right);

	Vector2();
	Vector2(const Vector2& vec);
	Vector2(const double& coor0, const double& coor1);
	~Vector2();
};