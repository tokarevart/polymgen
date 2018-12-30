#pragma once
#include "Vec3.h"
#include "Segment3.h"
#include "Ray3.h"

struct Vec3;
typedef Vec3 Point3;
struct Segment3;
typedef Segment3 Line3;
struct Ray3;

struct Track3
{
	Point3 origin;
	Vec3 velocity;

	//Point3& getOrigin() { return *origin; }

	Track3& operator= (const Track3& right);
	Track3  operator+ (const Vec3&   right) const;
	Track3  operator- (const Vec3&   right) const;
	Track3& operator+=(const Vec3&   right);
	Track3& operator-=(const Vec3&   right);

	explicit operator Line3() const;
	operator Ray3()  const;

	Track3();
	Track3(const Track3& ray);
	Track3(const Point3& origin, const Vec3& velocity);
	~Track3();
};