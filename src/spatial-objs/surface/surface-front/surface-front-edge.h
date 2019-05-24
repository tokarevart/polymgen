// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/surface/surface-face.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace surface {
namespace front {

class Edge
{
public:
    pmg::Edge* edge;
    vec3 normal;

    vec3 computeNormal();
    vec3 computeCenter();

    Edge( const surface::Face* relatedSurfaceFace, const pmg::Edge* edge );


private:
    // TODO: it's not good to store it
    surface::Face* m_relatedSurfaceFace;
};

} // namespace front
} // namespace surface
} // namespace pmg
