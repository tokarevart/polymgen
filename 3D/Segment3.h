#pragma once
#include "Vec3.h"

struct Segment3;
typedef Segment3 Line3;
struct Vec3;
typedef Vec3 Point3;

struct Segment3
{
	Point3* points;

	const double magnitude()    const;
	const double sqrMagnitude() const;
	
	Segment3& operator= (const Segment3& right);
	Segment3  operator+ (const Vec3&     right) const;
	Segment3  operator- (const Vec3&     right) const;
	Segment3& operator+=(const Vec3&     right);
	Segment3& operator-=(const Vec3&     right);
	
	Segment3();
	Segment3(const Segment3& segment);
	Segment3(const Point3& p0, const Point3& p1);
	~Segment3();
};