#include "Ray3.h"
#include "SpatialAlgorithms.h"


Ray3& Ray3::operator=(const Ray3& right)
{
	if (this == &right)
		return *this;

	origin = right.origin;
	direction = right.direction;
	return *this;
}

Ray3 Ray3::operator+(const Vec3& right) const
{
	return Ray3(origin + right, direction);
}

Ray3 Ray3::operator-(const Vec3& right) const
{
	return Ray3(origin - right, direction);
}

Ray3& Ray3::operator+=(const Vec3& right)
{
	origin += right;
	return *this;
}

Ray3& Ray3::operator-=(const Vec3& right)
{
	origin -= right;
	return *this;
}


Ray3::operator Line3() const
{
	return Line3(origin, origin + direction);
}

Ray3::operator Track3() const
{
	return Track3(origin, direction);
}


Ray3::Ray3() {}

Ray3::Ray3(const Ray3& ray)
{
	*this = ray;
}

Ray3::Ray3(const Point3& origin, const Vec3& direction)
	: origin(origin), direction(direction) {}

Ray3::~Ray3() {}