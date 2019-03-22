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




real_t vec3::dot(const vec3& vec0, const vec3& vec1)
{
    return vec0.coors[0] * vec1.coors[0] + vec0.coors[1] * vec1.coors[1] + vec0.coors[2] * vec1.coors[2];
}


vec3 vec3::cross(const vec3& vec0, const vec3& vec1)
{
    return vec3(vec0.coors[1] * vec1.coors[2] - vec0.coors[2] * vec1.coors[1],
               vec0.coors[2] * vec1.coors[0] - vec0.coors[0] * vec1.coors[2],
               vec0.coors[0] * vec1.coors[1] - vec0.coors[1] * vec1.coors[0]);
}


real_t vec3::mixed(const vec3 &vec0, const vec3 &vec1, const vec3 &vec2)
{
    return dot(cross(vec0, vec1), vec2);
}




real_t vec3::cos(const vec3& vec0, const vec3& vec1)
{
    if constexpr (std::is_same<real_t, float>())
    {
        return vec3::dot(vec0, vec1) / sqrtf(vec0.sqrMagnitude() * vec1.sqrMagnitude());
    }
    else
    {
        return vec3::dot(vec0, vec1) / sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
    }
}




real_t vec3::magnitude() const
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


real_t vec3::sqrMagnitude() const
{
    return dot(*this, *this);
}




vec3& vec3::normalize()
{
    real_t inv_magn = static_cast<real_t>(1.0) / magnitude();
    coors[0] *= inv_magn;
    coors[1] *= inv_magn;
    coors[2] *= inv_magn;
    return *this;
}


vec3& vec3::project(const vec3& vec)
{
    return *this = (dot(*this, vec) / vec.sqrMagnitude()) * vec;
}


vec3 vec3::project(const vec3& vec) const
{
    return vec3(*this).project(vec);
}


vec3& vec3::project(const vec3& plane_v0, const vec3& plane_v1)
{
    return *this -= vec3(*this).project(cross(plane_v0, plane_v1));
}


vec3 vec3::project(const vec3& plane_v0, const vec3& plane_v1) const
{
    return vec3(*this).project(plane_v0, plane_v1);
}




vec3& vec3::operator=(const vec3& right)
{
    coors[0] = right.coors[0];
    coors[1] = right.coors[1];
    coors[2] = right.coors[2];
    return *this;
}


vec3 vec3::operator+() const
{
    return vec3(*this);
}


vec3 vec3::operator-() const
{
    return vec3(-coors[0], -coors[1], -coors[2]);
}


vec3 vec3::operator+(const vec3& right) const
{
    return vec3(
        coors[0] + right.coors[0],
        coors[1] + right.coors[1],
        coors[2] + right.coors[2]);
}


vec3 vec3::operator-(const vec3& right) const
{
    return vec3(
        coors[0] - right.coors[0],
        coors[1] - right.coors[1],
        coors[2] - right.coors[2]);
}


vec3& vec3::operator+=(const vec3& right)
{
    coors[0] += right.coors[0];
    coors[1] += right.coors[1];
    coors[2] += right.coors[2];
    return *this;
}


vec3& vec3::operator-=(const vec3& right)
{
    coors[0] -= right.coors[0];
    coors[1] -= right.coors[1];
    coors[2] -= right.coors[2];
    return *this;
}


vec3& vec3::operator*=(real_t scalar)
{
    coors[0] *= scalar;
    coors[1] *= scalar;
    coors[2] *= scalar;
    return *this;
}


vec3& vec3::operator/=(real_t scalar)
{
    coors[0] /= scalar;
    coors[1] /= scalar;
    coors[2] /= scalar;
    return *this;
}


real_t& vec3::operator[](unsigned i)
{
    return coors[i];
}


const real_t& vec3::operator[](unsigned i) const
{
    return coors[i];
}


vec3 operator*(const vec3& vec, real_t scalar)
{
    return vec3(vec.coors[0] * scalar,
               vec.coors[1] * scalar,
               vec.coors[2] * scalar);
}


vec3 operator*(real_t scalar, const vec3& vec)
{
    return vec3(scalar * vec.coors[0],
               scalar * vec.coors[1],
               scalar * vec.coors[2]);
}


vec3 operator/(const vec3& vec, real_t scalar)
{
    real_t inv_scalar = static_cast<real_t>(1.0) / scalar;
    return vec3(vec.coors[0] * inv_scalar,
               vec.coors[1] * inv_scalar,
               vec.coors[2] * inv_scalar);
}




vec3::vec3()
{
    coors[0] = static_cast<real_t>(0.0);
    coors[1] = static_cast<real_t>(0.0);
    coors[2] = static_cast<real_t>(0.0);
}


vec3::vec3(const vec3& vec)
{
    *this = vec;
}


vec3::vec3(real_t coor0, real_t coor1, real_t coor2)
{
    coors[0] = coor0;
    coors[1] = coor1;
    coors[2] = coor2;
}
