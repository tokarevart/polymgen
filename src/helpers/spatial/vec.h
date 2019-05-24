// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <array>
#include "real-type.h"

// spt means spatial
namespace spt {

// TODO: try to make template and then make a benchmark
struct vec3
{
    using coordinate_t = real_t;

    std::array<coordinate_t, 3> x;

    static real_t dot  ( const vec3& vec0, const vec3& vec1 );
    static vec3   cross( const vec3& vec0, const vec3& vec1 );
    static real_t mixed( const vec3& vec0, const vec3& vec1, const vec3& vec2 );

    static real_t cos( const vec3& vec0, const vec3& vec1 );

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    vec3& normalize();
    vec3& project( const vec3& vec );
    vec3  project( const vec3& vec ) const;
    vec3& project( const vec3& plane_v0, const vec3& plane_v1 );
    vec3  project( const vec3& plane_v0, const vec3& plane_v1 ) const;

    vec3& operator=( const vec3& right );
    vec3  operator+() const;
    vec3  operator-() const;
    vec3  operator+( const vec3& right ) const;
    vec3  operator-( const vec3& right ) const;
    vec3& operator+=( const vec3& right );
    vec3& operator-=( const vec3& right );
    vec3& operator*=( real_t scalar );
    vec3& operator/=( real_t scalar );
    real_t&       operator[]( unsigned i );
    const real_t& operator[]( unsigned i ) const;

    vec3();
    vec3( const vec3& vec );
    vec3( coordinate_t x0, coordinate_t x1, coordinate_t x2 );
};

} // namespace spt

spt::vec3 operator*( const spt::vec3& vec, real_t scalar );
spt::vec3 operator*( real_t scalar, const spt::vec3& vec );
spt::vec3 operator/( const spt::vec3& vec, real_t scalar );
