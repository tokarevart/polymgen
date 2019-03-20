// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once

template <int Deg, typename Real = double>
constexpr Real cosDeg = static_cast<Real>(2.0);

#define SPEC_COSDEG(deg, value) \
    template<> \
    constexpr float  cosDeg<deg, float> = value ## f; \
    template<> \
    constexpr double cosDeg<deg, double> = value

// cosDeg specializations
SPEC_COSDEG(60,   0.5);
SPEC_COSDEG(70,   0.342020143325668733);
SPEC_COSDEG(80,   0.173648177666930348);
SPEC_COSDEG(90,   0.0);
SPEC_COSDEG(100, -0.173648177666930348);
SPEC_COSDEG(110, -0.342020143325668733);
SPEC_COSDEG(120, -0.5);
SPEC_COSDEG(130, -0.642787609686539326);
SPEC_COSDEG(140, -0.766044443118978035);
SPEC_COSDEG(150, -0.866025403784438646);
SPEC_COSDEG(160, -0.939692620785908384);
SPEC_COSDEG(170, -0.984807753012208059);
