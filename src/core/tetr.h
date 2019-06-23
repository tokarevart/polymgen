// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <algorithm>
#include "vert.h"
#include "../real-type.h"

#include "../definitions.h"


namespace pmg {

class Tetr {
    using FaceV = std::array<pmg::Vert*, 3>;

public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    // TODO: represent Tetr as a 4 Faces and then add method std::array<pmg::Vert*, 4> verts() const;
    std::array<pmg::Vert*, 4> verts;

    real_t computeVolume()  const;
    real_t computeQuality() const;

    bool contains(const Vert* vert) const;

    Tetr(const pmg::Vert* vert0,
         const pmg::Vert* vert1,
         const pmg::Vert* vert2,
         const pmg::Vert* vert3);
};

} // namespace pmg
