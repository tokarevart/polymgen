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
    using FaceV = std::array<Vert*, 3>;

public:
    Vert* verts[4];

    real_t computeVolume()  const;
    real_t computeQuality() const;

    bool contains( const Vert* vert ) const;

    Vert* findNot( FaceV face ) const;
    FaceV largestFace() const;

    static bool adjByFace( const Tetr* tetr0, const Tetr* tetr1 );

    Tetr( const Vert* vert0,
          const Vert* vert1,
          const Vert* vert2,
          const Vert* vert3 );
};

} // namespace pmg
