#pragma once

struct Vec3;
typedef Vec3 Point3;

struct Vec3
{
    double coors[3];

    static const double dotProduct  (const Vec3& vec0, const Vec3& vec1);
    static       Vec3   crossProduct(const Vec3& vec0, const Vec3& vec1);
    static const double mixedProduct(const Vec3& vec0, const Vec3& vec1, const Vec3& vec2);

    static const double cos(const Vec3& vec0, const Vec3& vec1);

    const double magnitude()    const;
    const double sqrMagnitude() const;

    Vec3& normalize();
    Vec3& project(const Vec3& vec);
    Vec3& project(const Vec3& plane_v0, const Vec3& plane_v1);

    double distanceToLine    (const Point3& p0, const Point3& p1)                   const;
    double distanceToSegment (const Point3& p0, const Point3& p1)                   const;
    double distanceToPlane   (const Point3& p0, const Point3& p1, const Point3& p2) const;
    double distanceToTriangle(const Point3& p0, const Point3& p1, const Point3& p2) const;
    
    Vec3& operator=(const Vec3& right);
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