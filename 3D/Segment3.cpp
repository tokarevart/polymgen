#include "Segment3.h"
#include "SpatialAlgorithms.h"


const double Segment3::magnitude() const
{
	return (points[1] - points[0]).magnitude();
}

const double Segment3::sqrMagnitude() const
{
	return (points[1] - points[0]).sqrMagnitude();
}


Segment3& Segment3::operator=(const Segment3& right)
{
	if (this == &right)
		return *this;

	points[0] = right.points[0];
	points[1] = right.points[1];
	return *this;
}

Segment3 Segment3::operator+(const Vec3& right) const
{
	return Segment3(points[0] + right, points[1] + right);
}

Segment3 Segment3::operator-(const Vec3& right) const
{
	return Segment3(points[0] - right, points[1] - right);
}

Segment3& Segment3::operator+=(const Vec3& right)
{
	points[0] += right;
	points[1] += right;
	return *this;
}

Segment3& Segment3::operator-=(const Vec3& right)
{
	points[0] -= right;
	points[1] -= right;
	return *this;
}


Segment3::Segment3() 
{
	points = new Point3[2];
}

Segment3::Segment3(const Segment3& segment)
{
	points = new Point3[2];
	*this = segment;
}

Segment3::Segment3(const Point3& p0, const Point3& p1)
{
	points = new Point3[2];
	points[0] = p0;
	points[1] = p1;
}

Segment3::~Segment3() 
{
	delete[] points;
}