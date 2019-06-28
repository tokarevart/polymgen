// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "../face.h"
#include "../../../helpers/spatial/vec.h"
#include "../../../real-type.h"

#include "../../../definitions.h"


namespace pmg::surface::front {

class Edge {
    using vec3 = spt::vec<3, real_t>;

public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    pmg::Edge* edge;
    vec3 normal;

    vec3 compute_normal();
    vec3 center();

    Edge(const surface::Face* relatedSurfaceFace, const pmg::Edge* edge);


private:
    // TODO: it's not good to store it
    surface::Face* m_related_surface_face;
};

} // namespace pmg::surface::front
