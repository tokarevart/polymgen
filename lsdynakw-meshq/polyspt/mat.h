#pragma once
#include "vec.h"


namespace spt {

template <std::size_t Dim, typename Real = typename spt::vec<Dim>::real_type>
struct mat;

template <typename Real>
struct mat<3, Real>
{
    static constexpr auto dim = 3;
    using real_type = Real;
    using line_type = spt::vec<3, real_type>;

    std::array<line_type, 3> x;

    line_type& operator[](std::uint8_t i)
    {
        return x[i];
    }
    const line_type& operator[](std::uint8_t i) const
    {
        return x[i];
    }

    mat() {}
    mat(const mat& other)
    {
        x = other.x;
    }
    mat(const std::array<line_type, 3>& x)
        : x(x) {}
    mat(const std::array<real_type, 9>& x)
    {
        this->x[0] = line_type(x[0], x[1], x[2]);
        this->x[1] = line_type(x[3], x[4], x[5]);
        this->x[2] = line_type(x[6], x[7], x[8]);
    }
    mat(const line_type& x0, const line_type& x1, const line_type& x2)
    {
        x[0] = x0;
        x[1] = x1;
        x[2] = x2;
    }
};

} // namespace spt
