// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "helpers/spatial/vec.h"

using spt::vec3;


real_t vec3::dot(const vec3& vec0, const vec3& vec1)
{
    return vec0.x[0] * vec1.x[0] + vec0.x[1] * vec1.x[1] + vec0.x[2] * vec1.x[2];
}


vec3 vec3::cross(const vec3& vec0, const vec3& vec1)
{
    return vec3(vec0.x[1] * vec1.x[2] - vec0.x[2] * vec1.x[1],
                vec0.x[2] * vec1.x[0] - vec0.x[0] * vec1.x[2],
                vec0.x[0] * vec1.x[1] - vec0.x[1] * vec1.x[0]);
}


real_t vec3::mixed(const vec3 &vec0, const vec3 &vec1, const vec3 &vec2)
{
    return dot(cross(vec0, vec1), vec2);
}




real_t vec3::cos(const vec3& vec0, const vec3& vec1)
{
    return vec3::dot(vec0, vec1) / std::sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
}




real_t vec3::magnitude() const
{
    return std::sqrt(sqrMagnitude());
}


real_t vec3::sqrMagnitude() const
{
    return dot(*this, *this);
}




vec3& vec3::normalize()
{
    real_t inv_magn = static_cast<real_t>(1.0) / magnitude();
    x[0] *= inv_magn;
    x[1] *= inv_magn;
    x[2] *= inv_magn;
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
    x[0] = right.x[0];
    x[1] = right.x[1];
    x[2] = right.x[2];
    return *this;
}


vec3 vec3::operator+() const
{
    return vec3(*this);
}


vec3 vec3::operator-() const
{
    return vec3(-x[0], -x[1], -x[2]);
}


vec3 vec3::operator+(const vec3& right) const
{
    return vec3(
        x[0] + right.x[0],
        x[1] + right.x[1],
        x[2] + right.x[2]);
}


vec3 vec3::operator-(const vec3& right) const
{
    return vec3(
        x[0] - right.x[0],
        x[1] - right.x[1],
        x[2] - right.x[2]);
}


vec3& vec3::operator+=(const vec3& right)
{
    x[0] += right.x[0];
    x[1] += right.x[1];
    x[2] += right.x[2];
    return *this;
}


vec3& vec3::operator-=(const vec3& right)
{
    x[0] -= right.x[0];
    x[1] -= right.x[1];
    x[2] -= right.x[2];
    return *this;
}


vec3& vec3::operator*=(real_t scalar)
{
    x[0] *= scalar;
    x[1] *= scalar;
    x[2] *= scalar;
    return *this;
}


vec3& vec3::operator/=(real_t scalar)
{
    x[0] /= scalar;
    x[1] /= scalar;
    x[2] /= scalar;
    return *this;
}


      real_t& vec3::operator[](unsigned i)       { return x[i]; }
const real_t& vec3::operator[](unsigned i) const { return x[i]; }


vec3 operator*(const vec3& vec, real_t scalar)
{
    return vec3(vec.x[0] * scalar,
                vec.x[1] * scalar,
                vec.x[2] * scalar);
}


vec3 operator*(real_t scalar, const vec3& vec)
{
    return vec3(scalar * vec.x[0],
                scalar * vec.x[1],
                scalar * vec.x[2]);
}


vec3 operator/(const vec3& vec, real_t scalar)
{
    real_t inv_scalar = static_cast<real_t>(1.0) / scalar;
    return vec3(vec.x[0] * inv_scalar,
                vec.x[1] * inv_scalar,
                vec.x[2] * inv_scalar);
}




vec3::vec3()
{
    x[0] = static_cast<real_t>(0.0);
    x[1] = static_cast<real_t>(0.0);
    x[2] = static_cast<real_t>(0.0);
}


vec3::vec3(const vec3& vec)
{
    *this = vec;
}


vec3::vec3(coordinate_t x0, coordinate_t x1, coordinate_t x2)
{
    x[0] = x0;
    x[1] = x1;
    x[2] = x2;
}
