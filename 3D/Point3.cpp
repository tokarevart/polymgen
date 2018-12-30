//#include "Vec3.h"
//#include "SpatialAlgorithms.h"
//
//#define DET(a, b, c, d) \
//		((a) * (d) - (b) * (c))
//
//#define BETWEEN(p0_coor, p1_coor, p) \
//		(std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)
//
//#define INSIDE_RECTANGLE(corner0, corner1, point) \
//		(BETWEEN(corner0.coors[0], corner1.coors[0], point.coors[0]) && \
//		 BETWEEN(corner0.coors[1], corner1.coors[1], point.coors[1]) && \
//		 BETWEEN(corner0.coors[2], corner1.coors[2], point.coors[2]))
//
//#define PI         3.141592653589793
//#define PI_DIV_180 0.017453292519943295
//
//
//double Vec3::distanceToLine(const Vec3& p0, const Vec3& p1) const
//{
//	return (spatialalgs::project(*this, { p0, p1 }) - *this).magnitude();
//}
//
//double Vec3::distanceToSegment(const Vec3& p0, const Vec3& p1) const
//{
//	Vec3 proj = spatialalgs::project(*this, { p0, p1 });
//	double sqr_magns[2];
//
//	const double EPS = (p1 - p0).magnitude() * 1e-6;
//	if (INSIDE_RECTANGLE(p0, p1, proj))
//		return (proj - *this).magnitude();
//
//	else if ((sqr_magns[0] = (p0 - *this).sqrMagnitude()) < (sqr_magns[1] = (p1 - *this).sqrMagnitude()))
//		return sqrt(sqr_magns[0]);
//
//	else
//		return sqrt(sqr_magns[1]);
//}
//
//
//Vec3& Vec3::operator=(const Vec3& right)
//{
//	if (this == &right)
//		return *this;
//
//	coors[0] = right.coors[0];
//	coors[1] = right.coors[1];
//	coors[2] = right.coors[2];
//	return *this;
//}
//
//Vec3 Vec3::operator+(const Vec3& right) const
//{
//	return Vec3(
//		coors[0] + right.coors[0],
//		coors[1] + right.coors[1],
//		coors[2] + right.coors[2]);
//}
//
//Vec3 Vec3::operator-(const Vec3& right) const
//{
//	return Vec3(
//		coors[0] - right.coors[0], 
//		coors[1] - right.coors[1], 
//		coors[2] - right.coors[2]);
//}
//
//Vec3 Vec3::operator-(const Vec3& right) const
//{
//	return Vec3(
//		coors[0] - right.coors[0], 
//		coors[1] - right.coors[1],
//		coors[2] - right.coors[2]);
//}
//
//Vec3& Vec3::operator+=(const Vec3& right)
//{
//	coors[0] += right.coors[0];
//	coors[1] += right.coors[1];
//	coors[2] += right.coors[2];
//	return *this;
//}
//
//Vec3& Vec3::operator-=(const Vec3& right)
//{
//	coors[0] -= right.coors[0];
//	coors[1] -= right.coors[1];
//	coors[2] -= right.coors[2];
//	return *this;
//}
//
//
//Vec3::operator Vec3() const
//{
//	return Vec3(coors[0], coors[1], coors[2]);
//}
//
//
//Vec3::Vec3()
//{
//	coors[0] = 0.0;
//	coors[1] = 0.0;
//	coors[2] = 0.0;
//}
//
//Vec3::Vec3(const Vec3& point)
//{
//	*this = point;
//}
//
//Vec3::Vec3(const double coor0, const double coor1, const double coor2)
//{
//	coors[0] = coor0;
//	coors[1] = coor1;
//	coors[2] = coor2;
//}
//
//Vec3::~Vec3() {}