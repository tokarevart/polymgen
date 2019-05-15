#pragma once
#include <stddef.h>
#include "../real-type.h"

namespace pmg {
namespace genparams {

struct Shell
{
    real_t cMinDis = static_cast<real_t>(0.2);
    real_t cMaxDef = static_cast<real_t>(0.3);
    real_t cDef    = static_cast<real_t>(0.5);

    size_t nSmoothIters         = 20;
    size_t nDelaunaySmoothIters =  3;
};

struct Volume
{
    real_t cMinDis = static_cast<real_t>(0.2);
    real_t cMaxDef = static_cast<real_t>(0.3);
    real_t cDef    = static_cast<real_t>(0.4);

    size_t nSmoothIters = 20;
};

struct Polyhedron
{
    Shell  shell  = Shell();
    Volume volume = Volume();
};

} // namespace genparams
} // namespace pmg
