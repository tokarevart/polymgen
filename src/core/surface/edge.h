// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include <vector>
#include "vert.h"
#include "../../real-type.h"

#include "../../definitions.h"


namespace pmg::surface {

class Edge {
public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    std::array<surface::Vert*, 2> verts;

    const std::vector<pmg::Edge*>& inner_edges() const;
    const std::vector<pmg::Vert*>& inner_verts() const;

    void segmentize(real_t preferredLen);

    real_t magnitude()    const;
    real_t sqr_magnitude() const;

    bool contains(const surface::Vert* sVert) const;
    bool contains(const     pmg::Edge* edge)  const;
    bool contains(const     pmg::Vert* vert)  const;

    Edge(const surface::Vert* vert0, const surface::Vert* vert1);


private:
    std::vector<pmg::Edge*> m_inner_edges;
    std::vector<pmg::Vert*> m_inner_verts;
};

} // namespace pmg::surface
