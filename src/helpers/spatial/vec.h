// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <array>
#include "../../real-type.h"

// spt means spatial
namespace spt {

// TODO: try to make template and then make a benchmark
struct vec3
{
    using coordinate_t = real_t;

    std::array<coordinate_t, 3> x = { static_cast<real_t>(0.0), };


    static real_t dot( const vec3& vec0, const vec3& vec1 )
    {
        return vec0.x[0] * vec1.x[0] + vec0.x[1] * vec1.x[1] + vec0.x[2] * vec1.x[2];
    }
    static vec3   cross( const vec3& vec0, const vec3& vec1 )
    {
        return vec3(vec0.x[1] * vec1.x[2] - vec0.x[2] * vec1.x[1],
                    vec0.x[2] * vec1.x[0] - vec0.x[0] * vec1.x[2],
                    vec0.x[0] * vec1.x[1] - vec0.x[1] * vec1.x[0]);
    }
    static real_t mixed( const vec3& vec0, const vec3& vec1, const vec3& vec2 )
    {
        return dot(cross(vec0, vec1), vec2);
    }

    static real_t cos( const vec3& vec0, const vec3& vec1 )
    {
        return vec3::dot(vec0, vec1) / std::sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
    }

    real_t magnitude() const
    {
        return std::sqrt(sqrMagnitude());
    }
    real_t sqrMagnitude() const
    {
        return dot(*this, *this);
    }

    vec3& normalize()
    {
        real_t inv_magn = static_cast<real_t>(1.0) / magnitude();
        x[0] *= inv_magn;
        x[1] *= inv_magn;
        x[2] *= inv_magn;
        return *this;
    }
    vec3& project( const vec3& vec )
    {
        return *this = vec * (dot(*this, vec) / vec.sqrMagnitude());
    }
    vec3  project( const vec3& vec ) const
    {
        return vec3(*this).project(vec);
    }
    vec3& project( const vec3& plane_v0, const vec3& plane_v1 )
    {
        return *this -= vec3(*this).project(cross(plane_v0, plane_v1));
    }
    vec3  project( const vec3& plane_v0, const vec3& plane_v1 ) const
    {
        return vec3(*this).project(plane_v0, plane_v1);
    }

    vec3& operator=( const vec3& right )
    {
        x = right.x; return *this;
    }
    vec3  operator+() const
    {
        return vec3(*this);
    }
    vec3  operator-() const
    {
        return vec3(-x[0], -x[1], -x[2]);
    }
    vec3  operator+( const vec3& right ) const
    {
        return vec3(x[0] + right.x[0],
                    x[1] + right.x[1],
                    x[2] + right.x[2]);
    }
    vec3  operator-( const vec3& right ) const
    {
        return vec3(x[0] - right.x[0],
                    x[1] - right.x[1],
                    x[2] - right.x[2]);
    }
    vec3  operator*( real_t scalar ) const
    {
        return vec3(x[0] * scalar,
                    x[1] * scalar,
                    x[2] * scalar);
    }
    vec3  operator/( real_t scalar ) const
    {
        real_t inv_scalar = static_cast<real_t>(1.0) / scalar;
        return vec3(x[0] * inv_scalar,
                    x[1] * inv_scalar,
                    x[2] * inv_scalar);
    }
    vec3& operator+=( const vec3& right )
    {
        x[0] += right.x[0];
        x[1] += right.x[1];
        x[2] += right.x[2];
        return *this;
    }
    vec3& operator-=( const vec3& right )
    {
        x[0] -= right.x[0];
        x[1] -= right.x[1];
        x[2] -= right.x[2];
        return *this;
    }
    vec3& operator*=( real_t scalar )
    {
        x[0] *= scalar;
        x[1] *= scalar;
        x[2] *= scalar;
        return *this;
    }
    vec3& operator/=( real_t scalar )
    {
        x[0] /= scalar;
        x[1] /= scalar;
        x[2] /= scalar;
        return *this;
    }
    real_t& operator[]( unsigned i )
    {
        return x[i];
    }
    const real_t& operator[]( unsigned i ) const
    {
        return x[i];
    }

    vec3() : x({ static_cast<real_t>(0.0), }) {}
    vec3( const vec3& vec )
    {
        x = vec.x;
    }
    vec3( const std::array<coordinate_t, 3>& x ) : x(x) {}
    vec3( coordinate_t x0, coordinate_t x1, coordinate_t x2 )
    {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
    }
};

} // namespace spt
