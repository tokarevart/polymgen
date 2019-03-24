// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include "real-type.h"


struct vec3
{
    real_t coors[3];

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
          real_t& operator[]( unsigned i );
    const real_t& operator[]( unsigned i ) const;

    vec3();
    vec3( const vec3& vec );
    vec3( real_t coor0, real_t coor1, real_t coor2 );
};

vec3 operator*( const vec3& vec, real_t scalar );
vec3 operator*( real_t scalar, const vec3& vec );
vec3 operator/( const vec3& vec, real_t scalar );
