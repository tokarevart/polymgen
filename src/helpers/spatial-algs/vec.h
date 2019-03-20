// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include "real-type.h"


struct Vec
{
    real_t coors[3];

    static real_t dot  ( const Vec& vec0, const Vec& vec1 );
    static Vec    cross( const Vec& vec0, const Vec& vec1 );
    static real_t mixed( const Vec& vec0, const Vec& vec1, const Vec& vec2 );

    static real_t cos( const Vec& vec0, const Vec& vec1 );

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    Vec& normalize();
    Vec& project( const Vec& vec );
    Vec  project( const Vec& vec ) const;
    Vec& project( const Vec& plane_v0, const Vec& plane_v1 );
    Vec  project( const Vec& plane_v0, const Vec& plane_v1 ) const;

    Vec& operator=( const Vec& right );
    Vec  operator+() const;
    Vec  operator-() const;
    Vec  operator+( const Vec& right ) const;
    Vec  operator-( const Vec& right ) const;
    Vec& operator+=( const Vec& right );
    Vec& operator-=( const Vec& right );
    Vec& operator*=( real_t scalar );
    Vec& operator/=( real_t scalar );
          real_t& operator[]( unsigned i );
    const real_t& operator[]( unsigned i ) const;

    Vec();
    Vec( const Vec& vec );
    Vec( real_t coor0, real_t coor1, real_t coor2 );
};

Vec operator*( const Vec& vec, real_t scalar );
Vec operator*( real_t scalar, const Vec& vec );
Vec operator/( const Vec& vec, real_t scalar );
