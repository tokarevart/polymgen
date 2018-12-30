#pragma once
//#include "Vec3.h"
//
//struct Vec3;
//
//struct Vec3
//{
//	double coors[3];
//
//	double distanceToLine    (const Vec3& p0, const Vec3& p1)                   const;
//	double distanceToSegment (const Vec3& p0, const Vec3& p1)                   const;
//	double distanceToPlane   (const Vec3& p0, const Vec3& p1, const Vec3& p2) const;
//	double distanceToTriangle(const Vec3& p0, const Vec3& p1, const Vec3& p2) const;
//	
//	Vec3& operator= (const Vec3& right);
//	Vec3  operator+ (const Vec3&   right) const;
//	Vec3  operator- (const Vec3&   right) const;
//	Vec3    operator- (const Vec3& right) const;
//	Vec3& operator+=(const Vec3& right);
//	Vec3& operator-=(const Vec3& right);
//
//	explicit operator Vec3() const;
//
//	Vec3();
//	Vec3(const Vec3& point);
//	Vec3(const double coor0, const double coor1, const double coor2);
//	~Vec3();
//};