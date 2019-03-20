// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <vector>
#include <memory>
#include "spatial-objs/vert.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace shell {

class Vert
{
public:
    pmg::Vert* attachedVert = nullptr;

    const Vec& pos() const;

          real_t& operator[](unsigned axis);
    const real_t& operator[](unsigned axis)   const;
    Vec operator-(const shell::Vert& other) const;
    Vec operator-(const   pmg::Vert& other) const;

    Vert();
    Vert(real_t coor0, real_t coor1, real_t coor2);
    Vert(const Vec& position);


private:
    std::unique_ptr<Vec> m_pos;
};

} // namespace shell
} // namespace pmg
