// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "../face.h"
#include "../../../helpers/spatial/vec.h"
#include "../../../real-type.h"

#include "../../../definitions.h"


namespace pmg::surface::front {

class Edge
{
public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    pmg::Edge* edge;
    spt::vec3 normal;

    spt::vec3 computeNormal();
    spt::vec3 computeCenter();

    Edge( const surface::Face* relatedSurfaceFace, const pmg::Edge* edge );


private:
    // TODO: it's not good to store it
    surface::Face* m_relatedSurfaceFace;
};

} // namespace pmg::surface::front
