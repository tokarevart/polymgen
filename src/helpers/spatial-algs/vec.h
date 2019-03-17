// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#pragma once
#include <cmath>


namespace tva {

struct Vec;
typedef Vec Point;

Vec operator*( const Vec& vec, double scalar );
Vec operator*( double scalar, const Vec& vec );
Vec operator/( const Vec& vec, double scalar );

struct Vec
{
    double coors[3];

    static double dot  ( const Vec& vec0, const Vec& vec1 );
    static Vec    cross( const Vec& vec0, const Vec& vec1 );
    static double mixed( const Vec& vec0, const Vec& vec1, const Vec& vec2 );

    static double cos( const Vec& vec0, const Vec& vec1 );

    double magnitude()    const;
    double sqrMagnitude() const;

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
    Vec& operator*=( double scalar );
    Vec& operator/=( double scalar );
          double& operator[]( unsigned i );
    const double& operator[]( unsigned i ) const;

    Vec();
    Vec( const Vec& vec );
    Vec( double coor0, double coor1, double coor2 );
};

} // namespace tva
