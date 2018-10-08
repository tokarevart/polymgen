#pragma once
#include <cmath>


class Vector3
{
private:
	double _coors[3];

public:
	static const double DotProduct(
		const Vector3& vec0, 
		const Vector3& vec1);
	static Vector3 CrossProduct(
		const Vector3& vec0, 
		const Vector3& vec1);
	static const double Cos(
		const Vector3& vec0, 
		const Vector3& vec1);
	static Vector3 Project(
		const Vector3& point, 
		const Vector3& linePoint0, 
		const Vector3& linePoint1);
	static const bool Project(
		Vector3& out, 
		const Vector3& point, 
		const Vector3& segmPoint0, 
		const Vector3& segmPoint1);
	static const bool RayIntersectTriangle(
		const Vector3& origin, 
		const Vector3& dir, 
		const Vector3& tr_v0, 
		const Vector3& tr_v1, 
		const Vector3& tr_v2);
	static const bool LineSegmentIntersectTriangle(
		const Vector3& segm_v0,
		const Vector3& segm_v1,
		const Vector3& tr_v0,
		const Vector3& tr_v1,
		const Vector3& tr_v2);

	const double Magnitude() const;
	const double SqrMagnitude() const;

	Vector3& Normalize();
	Vector3& Project(
		const Vector3& vec);
	Vector3& Project(
		const Vector3& planeVec0, 
		const Vector3& planeVec1);

	Vector3& operator=(const Vector3& vec);
	double& operator[](const int& axisIndex);
	Vector3 operator+() const;
	Vector3 operator-() const;
	Vector3 operator+(const Vector3& right) const;
	Vector3 operator-(const Vector3& right) const;
	Vector3& operator+=(const Vector3& right);
	Vector3& operator-=(const Vector3& right);
	Vector3& operator*=(const double& scalar);
	friend const Vector3 operator*(
		const Vector3& vec, 
		const double& scalar);
	friend const Vector3 operator*(
		const double& scalar, 
		const Vector3& vec);
	friend const Vector3 operator/(
		const Vector3& vec, 
		const double& scalar);

	Vector3();
	Vector3(
		const Vector3& vec);
	Vector3(
		const double& coor0, 
		const double& coor1, 
		const double& coor2);
	~Vector3();
};