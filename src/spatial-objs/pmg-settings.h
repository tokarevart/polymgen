#pragma once
#include <stddef.h>
#include "../real-type.h"

namespace pmg {
namespace settings {

struct Generation
{
    //
};

struct Optimization
{
    size_t nSmoothIters = 20;
    real_t minQuality   = static_cast<real_t>(0.3);
};

} // namespace settings
} // namespace pmg
