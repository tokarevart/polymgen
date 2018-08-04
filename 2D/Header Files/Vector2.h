#pragma once
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
	const double GetCoordinate(const int& axisIndex);
	void Rotate(const double& angle, AngleUnit unit);

	Vector2& operator=(const Vector2& vec);

	Vector2();
	Vector2(const Vector2& vec);
	Vector2(const double& coor0, const double& coor1);
	~Vector2();
};