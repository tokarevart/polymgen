// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "helpers/spatial-algs/vec.h"


#define DET(a, b, c, d) \
        ((a) * (d) - (b) * (c))

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) - EPS < (p) && (p) < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
            (BETWEEN(corner0.coors[0], corner1.coors[0], point.coors[0]) && \
             BETWEEN(corner0.coors[1], corner1.coors[1], point.coors[1]) && \
             BETWEEN(corner0.coors[2], corner1.coors[2], point.coors[2]))




real_t Vec::dot(const Vec& vec0, const Vec& vec1)
{
    return vec0.coors[0] * vec1.coors[0] + vec0.coors[1] * vec1.coors[1] + vec0.coors[2] * vec1.coors[2];
}


Vec Vec::cross(const Vec& vec0, const Vec& vec1)
{
    return Vec(vec0.coors[1] * vec1.coors[2] - vec0.coors[2] * vec1.coors[1],
               vec0.coors[2] * vec1.coors[0] - vec0.coors[0] * vec1.coors[2],
               vec0.coors[0] * vec1.coors[1] - vec0.coors[1] * vec1.coors[0]);
}


real_t Vec::mixed(const Vec &vec0, const Vec &vec1, const Vec &vec2)
{
    return dot(cross(vec0, vec1), vec2);
}




real_t Vec::cos(const Vec& vec0, const Vec& vec1)
{
    if constexpr (std::is_same<real_t, float>())
    {
        return Vec::dot(vec0, vec1) / sqrtf(vec0.sqrMagnitude() * vec1.sqrMagnitude());
    }
    else
    {
        return Vec::dot(vec0, vec1) / sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
    }
}




real_t Vec::magnitude() const
{
    if constexpr (std::is_same<real_t, float>())
    {
        return sqrtf(sqrMagnitude());
    }
    else
    {
        return sqrt(sqrMagnitude());
    }
}


real_t Vec::sqrMagnitude() const
{
    return dot(*this, *this);
}




Vec& Vec::normalize()
{
    real_t inv_magn = static_cast<real_t>(1.0) / magnitude();
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


Vec Vec::operator+() const
{
    return Vec(*this);
}


Vec Vec::operator-() const
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


Vec& Vec::operator*=(real_t scalar)
{
    coors[0] *= scalar;
    coors[1] *= scalar;
    coors[2] *= scalar;
    return *this;
}


Vec& Vec::operator/=(real_t scalar)
{
    coors[0] /= scalar;
    coors[1] /= scalar;
    coors[2] /= scalar;
    return *this;
}


real_t& Vec::operator[](unsigned i)
{
    return coors[i];
}


const real_t& Vec::operator[](unsigned i) const
{
    return coors[i];
}


Vec operator*(const Vec& vec, real_t scalar)
{
    return Vec(vec.coors[0] * scalar,
               vec.coors[1] * scalar,
               vec.coors[2] * scalar);
}


Vec operator*(real_t scalar, const Vec& vec)
{
    return Vec(scalar * vec.coors[0],
               scalar * vec.coors[1],
               scalar * vec.coors[2]);
}


Vec operator/(const Vec& vec, real_t scalar)
{
    real_t inv_scalar = static_cast<real_t>(1.0) / scalar;
    return Vec(vec.coors[0] * inv_scalar,
               vec.coors[1] * inv_scalar,
               vec.coors[2] * inv_scalar);
}




Vec::Vec()
{
    coors[0] = static_cast<real_t>(0.0);
    coors[1] = static_cast<real_t>(0.0);
    coors[2] = static_cast<real_t>(0.0);
}


Vec::Vec(const Vec& vec)
{
    *this = vec;
}


Vec::Vec(real_t coor0, real_t coor1, real_t coor2)
{
    coors[0] = coor0;
    coors[1] = coor1;
    coors[2] = coor2;
}
