#include "Vec3.h"
#include "SpatialAlgs.h"

#define DET(a, b, c, d) \
        ((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
            (BETWEEN(corner0.coors[0], corner1.coors[0], point.coors[0]) && \
             BETWEEN(corner0.coors[1], corner1.coors[1], point.coors[1]) && \
             BETWEEN(corner0.coors[2], corner1.coors[2], point.coors[2]))

#define PI         3.141592653589793
#define PI_DIV_180 0.017453292519943295


const double Vec3::dotProduct(const Vec3& vec0, const Vec3& vec1)
{
    return vec0.coors[0] * vec1.coors[0] + vec0.coors[1] * vec1.coors[1] + vec0.coors[2] * vec1.coors[2];
}

Vec3 Vec3::crossProduct(const Vec3& vec0, const Vec3& vec1)
{
    return Vec3(
        vec0.coors[1] * vec1.coors[2] - vec0.coors[2] * vec1.coors[1],
        vec0.coors[2] * vec1.coors[0] - vec0.coors[0] * vec1.coors[2],
        vec0.coors[0] * vec1.coors[1] - vec0.coors[1] * vec1.coors[0]);
}

const double Vec3::mixedProduct(const Vec3 &vec0, const Vec3 &vec1, const Vec3 &vec2)
{
    return dotProduct(crossProduct(vec0, vec1), vec2);
}


const double Vec3::cos(const Vec3& vec0, const Vec3& vec1)
{
    return Vec3::dotProduct(vec0, vec1) / sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
}


const double Vec3::magnitude() const
{
    return sqrt(sqrMagnitude());
}

const double Vec3::sqrMagnitude() const
{
    return dotProduct(*this, *this);
}


Vec3& Vec3::normalize()
{
    double inv_magn = 1.0 / magnitude();
    coors[0] *= inv_magn;
    coors[1] *= inv_magn;
    coors[2] *= inv_magn;
    return *this;
}

Vec3& Vec3::project(const Vec3& vec)
{
    return *this = (dotProduct(*this, vec) / vec.sqrMagnitude()) * vec;;
}

Vec3& Vec3::project(const Vec3& plane_v0, const Vec3& plane_v1)
{
    return *this -= Vec3(*this).project(crossProduct(plane_v0, plane_v1));;
}


double Vec3::distanceToLine(const Point3& p0, const Point3& p1) const
{
    return (spatialalgs::project(*this, p0, p1) - *this).magnitude();
}

double Vec3::distanceToSegment(const Point3& p0, const Point3& p1) const
{
    Point3 proj = spatialalgs::project(*this, p0, p1);
    double sqr_magns[2];

    const double EPS = (p1 - p0).magnitude() * 1e-6;
    if (INSIDE_RECTANGLE(p0, p1, proj))
        return (proj - *this).magnitude();
    else if ((sqr_magns[0] = (p0 - *this).sqrMagnitude()) < (sqr_magns[1] = (p1 - *this).sqrMagnitude()))
        return sqrt(sqr_magns[0]);
    else
        return sqrt(sqr_magns[1]);
}


Vec3& Vec3::operator=(const Vec3& right)
{
    if (this == &right)
        return *this;

    coors[0] = right.coors[0];
    coors[1] = right.coors[1];
    coors[2] = right.coors[2];
    return *this;
}

Vec3 Vec3::operator+() const
{
    return Vec3(*this);
}

Vec3 Vec3::operator-() const
{
    return Vec3(-coors[0], -coors[1], -coors[2]);
}

Vec3 Vec3::operator+(const Vec3& right) const
{
    return Vec3(
        coors[0] + right.coors[0], 
        coors[1] + right.coors[1], 
        coors[2] + right.coors[2]);
}

Vec3 Vec3::operator-(const Vec3& right) const
{
    return Vec3(
        coors[0] - right.coors[0], 
        coors[1] - right.coors[1], 
        coors[2] - right.coors[2]);
}

Vec3& Vec3::operator+=(const Vec3& right)
{
    coors[0] += right.coors[0];
    coors[1] += right.coors[1];
    coors[2] += right.coors[2];
    return *this;
}

Vec3& Vec3::operator-=(const Vec3& right)
{
    coors[0] -= right.coors[0];
    coors[1] -= right.coors[1];
    coors[2] -= right.coors[2];
    return *this;
}

Vec3& Vec3::operator*=(const double& scalar)
{
    coors[0] *= scalar;
    coors[1] *= scalar;
    coors[2] *= scalar;
    return *this;
}

Vec3& Vec3::operator/=(const double& scalar)
{
    coors[0] /= scalar;
    coors[1] /= scalar;
    coors[2] /= scalar;
    return *this;
}

const Vec3 operator*(const Vec3& vec, const double scalar)
{
    return Vec3(
        vec.coors[0] * scalar, 
        vec.coors[1] * scalar, 
        vec.coors[2] * scalar);
}

const Vec3 operator*(const double scalar, const Vec3& vec)
{
    return Vec3(
        scalar * vec.coors[0], 
        scalar * vec.coors[1], 
        scalar * vec.coors[2]);
}

const Vec3 operator/(const Vec3& vec, const double scalar)
{
    double inv_scalar = 1.0 / scalar;
    return Vec3(
        vec.coors[0] * inv_scalar, 
        vec.coors[1] * inv_scalar, 
        vec.coors[2] * inv_scalar);
}


Vec3::Vec3()
{
    coors[0] = 0.0;
    coors[1] = 0.0;
    coors[2] = 0.0;
}

Vec3::Vec3(const Vec3& vec)
{
    *this = vec;
}

Vec3::Vec3(const double coor0, const double coor1, const double coor2)
{
    coors[0] = coor0;
    coors[1] = coor1;
    coors[2] = coor2;
}

Vec3::~Vec3() {}
