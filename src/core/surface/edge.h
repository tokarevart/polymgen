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

    const std::vector<pmg::Edge*>& innerEdges() const;
    const std::vector<pmg::Vert*>& innerVerts() const;

    void segmentize(real_t preferredLen);

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    bool contains(const surface::Vert* sVert) const;
    bool contains(const     pmg::Edge* edge)  const;
    bool contains(const     pmg::Vert* vert)  const;

    Edge(const surface::Vert* vert0, const surface::Vert* vert1);


private:
    std::vector<pmg::Edge*> m_innerEdges;
    std::vector<pmg::Vert*> m_innerVerts;
};

} // namespace pmg::surface
