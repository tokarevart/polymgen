#pragma once
#include <cmath>

class Vec3
{
	double _coors[3];

public:
	static const double dotProduct  (const Vec3& vec0,  const Vec3& vec1);
	static       Vec3   crossProduct(const Vec3& vec0,  const Vec3& vec1);
	static const double mixedProduct(const Vec3& vec0,  const Vec3& vec1, const Vec3& vec2);

	static const double cos(const Vec3& vec0,  const Vec3& vec1);

	static       Vec3   project(const Vec3& point, const Vec3& line_p0, const Vec3& line_p1);
	static const bool   project(Vec3& out, const Vec3& point, const Vec3& segm_p0, const Vec3& segm_p1);

	static const bool   rayIntersectTriangle    (const Vec3& origin,  const Vec3& dir,     const Vec3& tr_p0, const Vec3& tr_p1, const Vec3& tr_p2);
	static const bool   segmentIntersectTriangle(const Vec3& segm_p0, const Vec3& segm_p1, const Vec3& tr_p0, const Vec3& tr_p1, const Vec3& tr_p2);
	
	static const double linesDistance   (const Vec3& l0_p0, const Vec3& l0_p1, const Vec3& l1_p0, const Vec3& l1_p1);
	static const double segmentsDistance(const Vec3& s0_p0, const Vec3& s0_p1, const Vec3& s1_p0, const Vec3& s1_p1);

	// CPA - Closest Point of Approach.
	static const double cpaTime    (const Vec3& start0, const Vec3& velocity0, const Vec3& start1, const Vec3& velocity1);
	static const double cpaDistance(const Vec3& start0, const Vec3& velocity0, const Vec3& start1, const Vec3& velocity1);

	const double magnitude()    const;
	const double sqrMagnitude() const;

	Vec3& normalize();
	Vec3& project(const Vec3& vec);
	Vec3& project(const Vec3& plane_v0, const Vec3& plane_v1);

	double distanceToLine   (const Vec3& line_p0, const Vec3& line_p1)                         const;
	double distanceToSegment(const Vec3& segm_p0, const Vec3& segm_p1)                         const;
	double distanceToPlane  (const Vec3& plane_p0, const Vec3& plane_p1, const Vec3& plane_p2) const;

	Vec3& operator=(const Vec3& vec);
	double& operator[](const int axisIndex);
	Vec3 operator+() const;
	Vec3 operator-() const;
	Vec3 operator+(const Vec3& right) const;
	Vec3 operator-(const Vec3& right) const;
	Vec3& operator+=(const Vec3& right);
	Vec3& operator-=(const Vec3& right);
	Vec3& operator*=(const double& scalar);
	Vec3& operator/=(const double& scalar);
	friend const Vec3 operator*(const Vec3& vec, const double scalar);
	friend const Vec3 operator*(const double scalar, const Vec3& vec);
	friend const Vec3 operator/(const Vec3& vec, const double scalar);

	Vec3();
	Vec3(const Vec3& vec);
	Vec3(const double coor0, const double coor1, const double coor2);
	~Vec3();
};