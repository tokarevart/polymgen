#pragma once
#include <cmath>

class Vector3;
typedef Vector3 Point3;

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
		const Point3& point,
		const Point3& line_p0,
		const Point3& line_p1);
	static const bool Project(
		Point3& out,
		const Point3& point,
		const Point3& segm_p0,
		const Point3& segm_p1);
	static const bool RayIntersectTriangle(
		const Point3& origin,
		const Vector3& dir, 
		const Point3& tr_p0,
		const Point3& tr_p1,
		const Point3& tr_p2);
	static const bool SegmentIntersectTriangle(
		const Point3& segm_p0,
		const Point3& segm_p1,
		const Point3& tr_p0,
		const Point3& tr_p1,
		const Point3& tr_p2);
	static const double LinesDistance(
		const Point3& l0_p0,
		const Point3& l0_p1,
		const Point3& l1_p0,
		const Point3& l1_p1);
	static const double SegmentsDistance(
		const Point3& s0_p0,
		const Point3& s0_p1,
		const Point3& s1_p0,
		const Point3& s1_p1);

	const double Magnitude() const;
	const double SqrMagnitude() const;

	Vector3& Normalize();
	Vector3& Project(
		const Vector3& vec);
	Vector3& Project(
		const Vector3& plane_v0, 
		const Vector3& plane_v1);

	double DistanceToLine(
		const Point3& line_p0,
		const Point3& line_p1) const;
	double DistanceToSegment(
		const Point3& segm_p0,
		const Point3& segm_p1) const;
	double DistanceToPlane(
		const Point3& plane_p0,
		const Point3& plane_p1,
		const Point3& plane_p2) const;

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