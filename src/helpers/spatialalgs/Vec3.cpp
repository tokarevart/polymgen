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




double tva::Vec3::dotProduct(const tva::Vec3& vec0, const tva::Vec3& vec1)
{
    return vec0.coors[0] * vec1.coors[0] + vec0.coors[1] * vec1.coors[1] + vec0.coors[2] * vec1.coors[2];
}


tva::Vec3 tva::Vec3::crossProduct(const tva::Vec3& vec0, const tva::Vec3& vec1)
{
    return tva::Vec3(
        vec0.coors[1] * vec1.coors[2] - vec0.coors[2] * vec1.coors[1],
        vec0.coors[2] * vec1.coors[0] - vec0.coors[0] * vec1.coors[2],
        vec0.coors[0] * vec1.coors[1] - vec0.coors[1] * vec1.coors[0]);
}


double tva::Vec3::mixedProduct(const tva::Vec3 &vec0, const tva::Vec3 &vec1, const tva::Vec3 &vec2)
{
    return dotProduct(crossProduct(vec0, vec1), vec2);
}




double tva::Vec3::cos(const tva::Vec3& vec0, const tva::Vec3& vec1)
{
    return tva::Vec3::dotProduct(vec0, vec1) / sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
}




double tva::Vec3::magnitude() const
{
    return sqrt(sqrMagnitude());
}


double tva::Vec3::sqrMagnitude() const
{
    return dotProduct(*this, *this);
}




tva::Vec3& tva::Vec3::normalize()
{
    double inv_magn = 1.0 / magnitude();
    coors[0] *= inv_magn;
    coors[1] *= inv_magn;
    coors[2] *= inv_magn;
    return *this;
}


tva::Vec3& tva::Vec3::project(const tva::Vec3& vec)
{
    return *this = (dotProduct(*this, vec) / vec.sqrMagnitude()) * vec;;
}


tva::Vec3& tva::Vec3::project(const tva::Vec3& plane_v0, const tva::Vec3& plane_v1)
{
    return *this -= tva::Vec3(*this).project(crossProduct(plane_v0, plane_v1));;
}




double tva::Vec3::distanceToLine(const tva::Point3& p0, const tva::Point3& p1) const
{
    return (spatialalgs::project(*this, p0, p1) - *this).magnitude();
}


double tva::Vec3::distanceToSegment(const tva::Point3& p0, const tva::Point3& p1) const
{
    tva::Point3 proj = spatialalgs::project(*this, p0, p1);
    double sqr_magns[2];

    const double EPS = (p1 - p0).magnitude() * 1e-6;
    if (INSIDE_RECTANGLE(p0, p1, proj))
        return (proj - *this).magnitude();
    else if ((sqr_magns[0] = (p0 - *this).sqrMagnitude()) < (sqr_magns[1] = (p1 - *this).sqrMagnitude()))
        return sqrt(sqr_magns[0]);
    else
        return sqrt(sqr_magns[1]);
}




tva::Vec3& tva::Vec3::operator=(const tva::Vec3& right)
{
    if (this == &right)
        return *this;

    coors[0] = right.coors[0];
    coors[1] = right.coors[1];
    coors[2] = right.coors[2];
    return *this;
}


tva::Vec3 tva::Vec3::operator+() const
{
    return tva::Vec3(*this);
}


tva::Vec3 tva::Vec3::operator-() const
{
    return tva::Vec3(-coors[0], -coors[1], -coors[2]);
}


tva::Vec3 tva::Vec3::operator+(const tva::Vec3& right) const
{
    return tva::Vec3(
        coors[0] + right.coors[0], 
        coors[1] + right.coors[1], 
        coors[2] + right.coors[2]);
}


tva::Vec3 tva::Vec3::operator-(const tva::Vec3& right) const
{
    return tva::Vec3(
        coors[0] - right.coors[0], 
        coors[1] - right.coors[1], 
        coors[2] - right.coors[2]);
}


tva::Vec3& tva::Vec3::operator+=(const tva::Vec3& right)
{
    coors[0] += right.coors[0];
    coors[1] += right.coors[1];
    coors[2] += right.coors[2];
    return *this;
}


tva::Vec3& tva::Vec3::operator-=(const tva::Vec3& right)
{
    coors[0] -= right.coors[0];
    coors[1] -= right.coors[1];
    coors[2] -= right.coors[2];
    return *this;
}


tva::Vec3& tva::Vec3::operator*=(double scalar)
{
    coors[0] *= scalar;
    coors[1] *= scalar;
    coors[2] *= scalar;
    return *this;
}


tva::Vec3& tva::Vec3::operator/=(double scalar)
{
    coors[0] /= scalar;
    coors[1] /= scalar;
    coors[2] /= scalar;
    return *this;
}


tva::Vec3 tva::operator*(const tva::Vec3& vec, double scalar)
{
    return tva::Vec3(
        vec.coors[0] * scalar, 
        vec.coors[1] * scalar, 
        vec.coors[2] * scalar);
}


tva::Vec3 tva::operator*(double scalar, const tva::Vec3& vec)
{
    return tva::Vec3(
        scalar * vec.coors[0], 
        scalar * vec.coors[1], 
        scalar * vec.coors[2]);
}


tva::Vec3 tva::operator/(const tva::Vec3& vec, double scalar)
{
    double inv_scalar = 1.0 / scalar;
    return tva::Vec3(
        vec.coors[0] * inv_scalar, 
        vec.coors[1] * inv_scalar, 
        vec.coors[2] * inv_scalar);
}




tva::Vec3::Vec3()
{
    coors[0] = 0.0;
    coors[1] = 0.0;
    coors[2] = 0.0;
}


tva::Vec3::Vec3(const tva::Vec3& vec)
{
    *this = vec;
}


tva::Vec3::Vec3(double coor0, double coor1, double coor2)
{
    coors[0] = coor0;
    coors[1] = coor1;
    coors[2] = coor2;
}


tva::Vec3::~Vec3() {}
