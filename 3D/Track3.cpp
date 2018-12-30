#include "Track3.h"
#include "SpatialAlgorithms.h"


Track3& Track3::operator=(const Track3& right)
{
	if (this == &right)
		return *this;

	origin = right.origin;
	velocity = right.velocity;
	return *this;
}

Track3 Track3::operator+(const Vec3& right) const
{
	return Ray3(origin + right, velocity);
}

Track3 Track3::operator-(const Vec3& right) const
{
	return Ray3(origin - right, velocity);
}

Track3& Track3::operator+=(const Vec3& right)
{
	origin += right;
	return *this;
}

Track3& Track3::operator-=(const Vec3& right)
{
	origin -= right;
	return *this;
}


Track3::operator Line3() const
{
	return Line3(origin, origin + velocity);
}

Track3::operator Ray3() const
{
	return Ray3(origin, velocity);
}


Track3::Track3() {}

Track3::Track3(const Track3& track)
{
	*this = track;
}

Track3::Track3(const Point3& origin, const Vec3& velocity)
	: origin(origin), velocity(velocity) {}

Track3::~Track3() {}