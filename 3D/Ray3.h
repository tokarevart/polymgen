#pragma once
#include "Vec3.h"
#include "Segment3.h"
#include "Track3.h"

struct Vec3;
typedef Vec3 Point3;
struct Segment3;
typedef Segment3 Line3;
struct Track3;

struct Ray3
{
	Point3 origin;
	Vec3 direction;

	Ray3& operator= (const Ray3& right);
	Ray3  operator+ (const Vec3& right) const;
	Ray3  operator- (const Vec3& right) const;
	Ray3& operator+=(const Vec3& right);
	Ray3& operator-=(const Vec3& right);

	explicit operator Line3()  const;
	operator Track3() const;

	Ray3();
	Ray3(const Ray3& ray);
	Ray3(const Point3& origin, const Vec3& direction);
	~Ray3();
};