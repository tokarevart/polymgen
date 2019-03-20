#include "real-type.h"

real_t sqrtReal( real_t value )
{
    if constexpr (std::is_same<real_t, float>())
    {
        return sqrtf(value);
    }
    else
    {
        return sqrt(value);
    }
}

real_t acosReal( real_t value )
{
    if constexpr (std::is_same<real_t, float>())
    {
        return acosf(value);
    }
    else
    {
        return acos(value);
    }
}

real_t roundReal( real_t value )
{
    if constexpr (std::is_same<real_t, float>())
    {
        return roundf(value);
    }
    else
    {
        return round(value);
    }
}
