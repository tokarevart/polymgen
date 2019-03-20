// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <algorithm>
#include "spatial-objs/vert.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Tetr
{
public:
    Vert* verts[4];

    real_t computeVolume()  const;
    real_t computeQuality() const;

    Tetr( const Vert* vert0,
          const Vert* vert1,
          const Vert* vert2,
          const Vert* vert3 );
};

} // namespace pmg
