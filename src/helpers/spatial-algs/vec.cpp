// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#include "helpers/spatial-algs/vec.h"


using tva::Vec;
using tva::Point;


#define DET(a, b, c, d) \
        ((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
            (BETWEEN(corner0.coors[0], corner1.coors[0], point.coors[0]) && \
             BETWEEN(corner0.coors[1], corner1.coors[1], point.coors[1]) && \
             BETWEEN(corner0.coors[2], corner1.coors[2], point.coors[2]))




double Vec::dot(const Vec& vec0, const Vec& vec1)
{
    return vec0.coors[0] * vec1.coors[0] + vec0.coors[1] * vec1.coors[1] + vec0.coors[2] * vec1.coors[2];
}


Vec Vec::cross(const Vec& vec0, const Vec& vec1)
{
    return Vec(
        vec0.coors[1] * vec1.coors[2] - vec0.coors[2] * vec1.coors[1],
        vec0.coors[2] * vec1.coors[0] - vec0.coors[0] * vec1.coors[2],
        vec0.coors[0] * vec1.coors[1] - vec0.coors[1] * vec1.coors[0]);
}


double Vec::mixed(const Vec &vec0, const Vec &vec1, const Vec &vec2)
{
    return dot(cross(vec0, vec1), vec2);
}




double Vec::cos(const Vec& vec0, const Vec& vec1)
{
    return Vec::dot(vec0, vec1) / sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
}




double Vec::magnitude() const
{
    return sqrt(sqrMagnitude());
}


double Vec::sqrMagnitude() const
{
    return dot(*this, *this);
}




Vec& Vec::normalize()
{
    double inv_magn = 1.0 / magnitude();
    coors[0] *= inv_magn;
    coors[1] *= inv_magn;
    coors[2] *= inv_magn;
    return *this;
}


Vec& Vec::project(const Vec& vec)
{
    return *this = (dot(*this, vec) / vec.sqrMagnitude()) * vec;
}


Vec Vec::project(const Vec& vec) const
{
    return Vec(*this).project(vec);
}


Vec& Vec::project(const Vec& plane_v0, const Vec& plane_v1)
{
    return *this -= Vec(*this).project(cross(plane_v0, plane_v1));
}


Vec Vec::project(const Vec& plane_v0, const Vec& plane_v1) const
{
    return Vec(*this).project(plane_v0, plane_v1);
}




Vec& Vec::operator=(const Vec& right)
{
    coors[0] = right.coors[0];
    coors[1] = right.coors[1];
    coors[2] = right.coors[2];
    return *this;
}


tva::Vec Vec::operator+() const
{
    return Vec(*this);
}


tva::Vec Vec::operator-() const
{
    return Vec(-coors[0], -coors[1], -coors[2]);
}


Vec Vec::operator+(const Vec& right) const
{
    return Vec(
        coors[0] + right.coors[0],
        coors[1] + right.coors[1],
        coors[2] + right.coors[2]);
}


Vec Vec::operator-(const Vec& right) const
{
    return Vec(
        coors[0] - right.coors[0],
        coors[1] - right.coors[1],
        coors[2] - right.coors[2]);
}


Vec& Vec::operator+=(const Vec& right)
{
    coors[0] += right.coors[0];
    coors[1] += right.coors[1];
    coors[2] += right.coors[2];
    return *this;
}


Vec& Vec::operator-=(const Vec& right)
{
    coors[0] -= right.coors[0];
    coors[1] -= right.coors[1];
    coors[2] -= right.coors[2];
    return *this;
}


Vec& Vec::operator*=(double scalar)
{
    coors[0] *= scalar;
    coors[1] *= scalar;
    coors[2] *= scalar;
    return *this;
}


Vec& Vec::operator/=(double scalar)
{
    coors[0] /= scalar;
    coors[1] /= scalar;
    coors[2] /= scalar;
    return *this;
}


double& Vec::operator[](unsigned i)
{
    return coors[i];
}


const double& Vec::operator[](unsigned i) const
{
    return coors[i];
}


Vec tva::operator*(const Vec& vec, double scalar)
{
    return Vec(
        vec.coors[0] * scalar,
        vec.coors[1] * scalar,
        vec.coors[2] * scalar);
}


Vec tva::operator*(double scalar, const Vec& vec)
{
    return Vec(
        scalar * vec.coors[0],
        scalar * vec.coors[1],
        scalar * vec.coors[2]);
}


Vec tva::operator/(const Vec& vec, double scalar)
{
    double inv_scalar = 1.0 / scalar;
    return Vec(
        vec.coors[0] * inv_scalar,
        vec.coors[1] * inv_scalar,
        vec.coors[2] * inv_scalar);
}




Vec::Vec()
{
    coors[0] = 0.0;
    coors[1] = 0.0;
    coors[2] = 0.0;
}


Vec::Vec(const Vec& vec)
{
    *this = vec;
}


Vec::Vec(double coor0, double coor1, double coor2)
{
    coors[0] = coor0;
    coors[1] = coor1;
    coors[2] = coor2;
}
