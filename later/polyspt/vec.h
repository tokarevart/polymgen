// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <numeric>
#include <array>

// spt means spatial
namespace spt {

template <std::size_t Dim, typename Real = double>
struct vec
{
    static constexpr auto dim = Dim;
    using real_type = Real;

    std::array<real_type, Dim> x = { static_cast<real_type>(0.0), };


    static real_type dot(const vec& vec0, const vec& vec1)
    {
        return std::inner_product(vec0.x.begin(), vec0.x.end(), vec1.x.begin(), static_cast<real_type>(0.0));
    }
   
    static real_type cos(const vec& vec0, const vec& vec1)
    {
        return vec::dot(vec0, vec1) / std::sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
    }

    real_type magnitude() const
    {
        return std::sqrt(sqrMagnitude());
    }
    real_type sqrMagnitude() const
    {
        return dot(*this, *this);
    }

    vec& normalize()
    {
        real_type inv_magn = static_cast<real_type>(1.0) / magnitude();
        for (auto& coor : x) coor *= inv_magn;
        return *this;
    }
    vec& project( const vec& vec )
    {
        return *this = vec * (dot(*this, vec) / vec.sqrMagnitude());
    }
    vec project( const vec& vec ) const
    {
        return vec(*this).project(vec);
    }

    vec& operator=(const vec& right)
    {
        x = right.x; 
        return *this;
    }
    vec operator-() const
    {
        vec res;
        for (std::size_t i = 0; i < Dim; i++)
            res.x[i] = -x[i];
        return res;
    }
    vec operator+(const vec& right) const
    {
        vec res;
        for (std::size_t i = 0; i < Dim; i++)
            res.x[i] = x[i] + right.x[i];
        return res;
    }
    vec operator-(const vec& right) const
    {
        vec res;
        for (std::size_t i = 0; i < Dim; i++)
            res.x[i] = x[i] - right.x[i];
        return res;
    }
    vec operator*(real_type scalar) const
    {
        vec res;
        for (std::size_t i = 0; i < Dim; i++)
            res.x[i] = x[i] * scalar;
        return res;
    }
    vec operator/(real_type scalar) const
    {
        vec res;
        for (std::size_t i = 0; i < Dim; i++)
            res.x[i] = x[i] / scalar;
        return res;
    }
    vec& operator+=(const vec& right)
    {
        for (std::size_t i = 0; i < Dim; i++)
            x[i] += right.x[i];
        return *this;
    }
    vec& operator-=(const vec& right)
    {
        for (std::size_t i = 0; i < Dim; i++)
            x[i] -= right.x[i];
        return *this;
    }
    vec& operator*=(real_type scalar)
    {
        for (std::size_t i = 0; i < Dim; i++)
            x[i] *= scalar;
        return *this;
    }
    vec& operator/=(real_type scalar)
    {
        for (std::size_t i = 0; i < Dim; i++)
            x[i] /= scalar;
        return *this;
    }
    real_type& operator[](std::size_t i)
    {
        return x[i];
    }
    const real_type& operator[](std::size_t i) const
    {
        return x[i];
    }

    vec() 
        : x({ static_cast<real_type>(0.0), }) {}
    vec(const vec& vec)
    {
        x = vec.x;
    }
    vec(const std::array<real_type, Dim>& x) 
        : x(x) {}
    template <typename... Reals>
    vec(Reals... xs) 
        : x({ xs... }) {}
};


template <typename Real>
struct vec<3, Real>
{
    using real_type = Real;

    std::array<real_type, 3> x = { static_cast<real_type>(0.0), };


    static real_type dot(const vec& vec0, const vec& vec1)
    {
        return vec0.x[0] * vec1.x[0] + vec0.x[1] * vec1.x[1] + vec0.x[2] * vec1.x[2];
    }
    static vec cross(const vec& vec0, const vec& vec1)
    {
        return vec(vec0.x[1] * vec1.x[2] - vec0.x[2] * vec1.x[1],
            vec0.x[2] * vec1.x[0] - vec0.x[0] * vec1.x[2],
            vec0.x[0] * vec1.x[1] - vec0.x[1] * vec1.x[0]);
    }
    static real_type mixed(const vec& vec0, const vec& vec1, const vec& vec2)
    {
        return dot(cross(vec0, vec1), vec2);
    }

    static real_type cos(const vec& vec0, const vec& vec1)
    {
        return vec::dot(vec0, vec1) / std::sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
    }

    real_type magnitude() const
    {
        return std::sqrt(sqrMagnitude());
    }
    real_type sqrMagnitude() const
    {
        return dot(*this, *this);
    }

    vec& normalize()
    {
        real_type inv_magn = static_cast<real_type>(1.0) / magnitude();
        x[0] *= inv_magn;
        x[1] *= inv_magn;
        x[2] *= inv_magn;
        return *this;
    }
    vec& project(const vec& vec)
    {
        return *this = vec * (dot(*this, vec) / vec.sqrMagnitude());
    }
    vec project(const vec& vec) const
    {
        return vec(*this).project(vec);
    }
    vec& project(const vec& plane_v0, const vec& plane_v1)
    {
        return *this -= vec(*this).project(cross(plane_v0, plane_v1));
    }
    vec project(const vec& plane_v0, const vec& plane_v1) const
    {
        return vec(*this).project(plane_v0, plane_v1);
    }

    vec& operator=(const vec& right)
    {
        x = right.x;
        return *this;
    }
    vec operator-() const
    {
        return vec(-x[0], -x[1], -x[2]);
    }
    vec operator+(const vec& right) const
    {
        return vec(x[0] + right.x[0],
                   x[1] + right.x[1],
                   x[2] + right.x[2]);
    }
    vec operator-(const vec& right) const
    {
        return vec(x[0] - right.x[0],
                   x[1] - right.x[1],
                   x[2] - right.x[2]);
    }
    vec operator*(real_type scalar) const
    {
        return vec(x[0] * scalar,
                   x[1] * scalar,
                   x[2] * scalar);
    }
    vec operator/(real_type scalar) const
    {
        return vec(x[0] / scalar,
                   x[1] / scalar,
                   x[2] / scalar);
    }
    vec& operator+=(const vec& right)
    {
        x[0] += right.x[0];
        x[1] += right.x[1];
        x[2] += right.x[2];
        return *this;
    }
    vec& operator-=(const vec& right)
    {
        x[0] -= right.x[0];
        x[1] -= right.x[1];
        x[2] -= right.x[2];
        return *this;
    }
    vec& operator*=(real_type scalar)
    {
        x[0] *= scalar;
        x[1] *= scalar;
        x[2] *= scalar;
        return *this;
    }
    vec& operator/=(real_type scalar)
    {
        x[0] /= scalar;
        x[1] /= scalar;
        x[2] /= scalar;
        return *this;
    }
    real_type& operator[](std::size_t i)
    {
        return x[i];
    }
    const real_type& operator[](std::size_t i) const
    {
        return x[i];
    }

    vec()
        : x({ static_cast<real_type>(0.0), }) {}
    vec(const vec& vec)
    {
        x = vec.x;
    }
    vec(const std::array<real_type, 3>& x)
        : x(x) {}
    vec(real_type x0, real_type x1, real_type x2)
    {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
    }
};


template <typename Real>
struct vec<2, Real>
{
    using real_type = Real;

    std::array<real_type, 2> x = { static_cast<real_type>(0.0), };


    static real_type dot(const vec& vec0, const vec& vec1)
    {
        return vec0.x[0] * vec1.x[0] + vec0.x[1] * vec1.x[1];
    }
    static real_type cross(const vec& vec0, const vec& vec1)
    {
        return vec0[0] * vec1[1] - vec0[1] * vec1[0];
    }

    static real_type cos(const vec& vec0, const vec& vec1)
    {
        return vec::dot(vec0, vec1) / std::sqrt(vec0.sqrMagnitude() * vec1.sqrMagnitude());
    }

    real_type magnitude() const
    {
        return std::sqrt(sqrMagnitude());
    }
    real_type sqrMagnitude() const
    {
        return dot(*this, *this);
    }

    vec& normalize()
    {
        real_type inv_magn = static_cast<real_type>(1.0) / magnitude();
        x[0] *= inv_magn;
        x[1] *= inv_magn;
        return *this;
    }
    vec& project(const vec& vec)
    {
        return *this = vec * (dot(*this, vec) / vec.sqrMagnitude());
    }
    vec project(const vec& vec) const
    {
        return vec(*this).project(vec);
    }

    vec& operator=(const vec& right)
    {
        x = right.x;
        return *this;
    }
    vec operator-() const
    {
        return vec(-x[0], -x[1]);
    }
    vec operator+(const vec& right) const
    {
        return vec(x[0] + right.x[0],
                   x[1] + right.x[1]);
    }
    vec operator-(const vec& right) const
    {
        return vec(x[0] - right.x[0],
                   x[1] - right.x[1]);
    }
    vec operator*(real_type scalar) const
    {
        return vec(x[0] * scalar,
                   x[1] * scalar);
    }
    vec operator/(real_type scalar) const
    {
        real_type inv_scalar = static_cast<real_type>(1.0) / scalar;
        return vec(x[0] * inv_scalar,
                   x[1] * inv_scalar);
    }
    vec& operator+=(const vec& right)
    {
        x[0] += right.x[0];
        x[1] += right.x[1];
        return *this;
    }
    vec& operator-=(const vec& right)
    {
        x[0] -= right.x[0];
        x[1] -= right.x[1];
        return *this;
    }
    vec& operator*=(real_type scalar)
    {
        x[0] *= scalar;
        x[1] *= scalar;
        return *this;
    }
    vec& operator/=(real_type scalar)
    {
        x[0] /= scalar;
        x[1] /= scalar;
        return *this;
    }
    real_type& operator[](std::size_t i)
    {
        return x[i];
    }
    const real_type& operator[](std::size_t i) const
    {
        return x[i];
    }

    vec()
        : x({ static_cast<real_type>(0.0), }) {}
    vec(const vec& vec)
    {
        x = vec.x;
    }
    vec(const std::array<real_type, 2>& x)
        : x(x) {}
    vec(real_type x0, real_type x1)
    {
        x[0] = x0;
        x[1] = x1;
    }
};

} // namespace spt
