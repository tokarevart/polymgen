#pragma once
#include <cstddef>
#include "../../real-type.h"
#include "../vert.h"


namespace pmg::front {

class Vert {
public:
    pmg::Vert* x;

    Vert(const pmg::Vert* vert) : x{ const_cast<pmg::Vert*>(vert) } {}
};

} // namespace pmg::front
