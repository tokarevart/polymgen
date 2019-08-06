// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cmath>
#include <cstdint>
#include <array>

namespace spt {

template <std::size_t Dim, typename ValueType = double>
struct vec;

template <typename Real> using vec2 = vec<2, Real>;
using vec2f = vec2<float>;
using vec2d = vec2<double>;
template <typename Real> using vec3 = vec<3, Real>;
using vec3f = vec3<float>;
using vec3d = vec3<double>;


template <typename ValueType>
struct vec<3, ValueType> {
    static constexpr std::size_t dim = 3;
    using value_type = ValueType;

    std::array<value_type, 3> x{
        static_cast<value_type>(0),
        static_cast<value_type>(0),
        static_cast<value_type>(0)
    };

    value_type magnitude() const {
        return std::sqrt(sqr_magnitude());
    }
    value_type sqr_magnitude() const {
        return dot(*this, *this);
    }

    vec& normalize() {
        auto inv_magn = static_cast<value_type>(1) / magnitude();
        x[0] *= inv_magn;
        x[1] *= inv_magn;
        x[2] *= inv_magn;
        return *this;
    }
    vec& project(const vec& v) {
        return *this = v * (dot(*this, v) / v.sqr_magnitude());
    }
    vec project(const vec& v) const {
        return vec(*this).project(v);
    }
    vec& project(const vec& plane_v0, const vec& plane_v1) {
        return *this -= vec(*this).project(cross(plane_v0, plane_v1));
    }
    vec project(const vec& plane_v0, const vec& plane_v1) const {
        return vec(*this).project(plane_v0, plane_v1);
    }

    vec& operator=(const vec& right) {
        x = right.x;
        return *this;
    }
    vec operator-() const {
        return { -x[0], -x[1], -x[2] };
    }
    vec operator+(const vec& right) const {
        return {
            x[0] + right.x[0],
            x[1] + right.x[1],
            x[2] + right.x[2] };
    }
    vec operator-(const vec& right) const {
        return {
            x[0] - right.x[0],
            x[1] - right.x[1],
            x[2] - right.x[2] };
    }
    vec operator*(value_type scalar) const {
        return {
            x[0] * scalar,
            x[1] * scalar,
            x[2] * scalar };
    }
    vec operator/(value_type scalar) const {
        return {
            x[0] / scalar,
            x[1] / scalar,
            x[2] / scalar };
    }
    vec& operator+=(const vec& right) {
        x[0] += right.x[0];
        x[1] += right.x[1];
        x[2] += right.x[2];
        return *this;
    }
    vec& operator-=(const vec& right) {
        x[0] -= right.x[0];
        x[1] -= right.x[1];
        x[2] -= right.x[2];
        return *this;
    }
    vec& operator*=(value_type scalar) {
        x[0] *= scalar;
        x[1] *= scalar;
        x[2] *= scalar;
        return *this;
    }
    vec& operator/=(value_type scalar) {
        x[0] /= scalar;
        x[1] /= scalar;
        x[2] /= scalar;
        return *this;
    }
    value_type& operator[](std::uint8_t i) {
        return x[i];
    }
    const value_type& operator[](std::uint8_t i) const {
        return x[i];
    }

    vec() {}
    vec(const vec& other) {
        x = other.x;
    }
    vec(const std::array<value_type, 3>& x)
        : x{ x } {}
    vec(value_type x0, value_type x1, value_type x2)
        : x{ x0, x1, x2 } {}
};


template <typename Real>
struct vec<2, Real> {
    static constexpr std::size_t dim = 2;
    using value_type = Real;

    std::array<value_type, 2> x{
        static_cast<value_type>(0),
        static_cast<value_type>(0)
    };

    value_type magnitude() const {
        return std::sqrt(sqr_magnitude());
    }
    value_type sqr_magnitude() const {
        return dot(*this, *this);
    }

    vec& normalize() {
        auto inv_magn = static_cast<value_type>(1.0) / magnitude();
        x[0] *= inv_magn;
        x[1] *= inv_magn;
        return *this;
    }
    vec& project(const vec& v) {
        return *this = v * (dot(*this, v) / v.sqr_magnitude());
    }
    vec project(const vec& v) const {
        return vec(*this).project(vec);
    }

    vec& operator=(const vec& right) {
        x = right.x;
        return *this;
    }
    vec operator-() const {
        return { -x[0], -x[1] };
    }
    vec operator+(const vec& right) const {
        return { x[0] + right.x[0], x[1] + right.x[1] };
    }
    vec operator-(const vec& right) const {
        return { x[0] - right.x[0], x[1] - right.x[1] };
    }
    vec operator*(value_type scalar) const {
        return { x[0] * scalar, x[1] * scalar };
    }
    vec operator/(value_type scalar) const {
        return { x[0] / scalar, x[1] / scalar };
    }
    vec& operator+=(const vec& right) {
        x[0] += right.x[0];
        x[1] += right.x[1];
        return *this;
    }
    vec& operator-=(const vec& right) {
        x[0] -= right.x[0];
        x[1] -= right.x[1];
        return *this;
    }
    vec& operator*=(value_type scalar) {
        x[0] *= scalar;
        x[1] *= scalar;
        return *this;
    }
    vec& operator/=(value_type scalar) {
        x[0] /= scalar;
        x[1] /= scalar;
        return *this;
    }
    value_type& operator[](std::uint8_t i) {
        return x[i];
    }
    const value_type& operator[](std::uint8_t i) const {
        return x[i];
    }

    vec() {}
    vec(const vec& other) {
        x = other.x;
    }
    vec(const std::array<value_type, 2>& x)
        : x{ x } {}
    vec(value_type x0, value_type x1)
        : x{ x0, x1 } {}
};

template <std::size_t N, typename ValueType>
vec<N, ValueType> operator*(ValueType scalar, const vec<N, ValueType>& v) {
    return v * scalar;
}

} // namespace spt
