//#include "Line3.h"
//
//
//Line3& Line3::operator=(const Line3& right)
//{
//	if (this == &right)
//		return *this;
//
//	points[0] = right.points[0];
//	points[1] = right.points[1];
//	return *this;
//}
//
//Line3 Line3::operator+(const Vec3& right) const
//{
//	return Line3(points[0] + right, points[1] + right);
//}
//
//Line3 Line3::operator-(const Vec3& right) const
//{
//	return Line3(points[0] - right, points[1] - right);
//}
//
//Line3& Line3::operator+=(const Vec3& right)
//{
//	points[0] += right;
//	points[1] += right;
//	return *this;
//}
//
//Line3& Line3::operator-=(const Vec3& right)
//{
//	points[0] -= right;
//	points[1] -= right;
//	return *this;
//}
//
//Line3::operator Segment3() const
//{
//	return Segment3(points[0], points[1]);
//}
//
//Line3::Line3() {}
//
//Line3::Line3(const Line3& line)
//{
//	*this = line;
//}
//
//Line3::Line3(const Vec3& p0, const Vec3& p1)
//{
//	points[0] = p0;
//	points[1] = p1;
//}
//
//Line3::~Line3() {}