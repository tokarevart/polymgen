// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <algorithm>
#include "spatial-objs/vertex.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Tetr
{
public:
    Vertex* verts[4];

    real_t computeVolume()  const;
    real_t computeQuality() const;

    Tetr( const Vertex* vert0,
          const Vertex* vert1,
          const Vertex* vert2,
          const Vertex* vert3 );
};

} // namespace pmg
