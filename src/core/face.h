// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include "edge.h"
#include "vert.h"
#include "../helpers/spatial/vec.h"
#include "../real-type.h"

#include "../definitions.h"


namespace pmg {

// TODO: generalize geom structures with templates (and no interfaces); use Polytope and so on; with static and dynamic num of subpolytopes
class Face {
    using vec3 = spt::vec<3, real_t>;

public:
    std::array<pmg::Edge*, 3> edges;

    // TODO: add method std::array<pmg::Vert*, 3> verts() const;

    vec3 center() const;
    real_t quality() const;
    real_t area() const;
    
    pmg::Vert* find_vert_not(const pmg::Edge* edge) const;
    pmg::Edge* find_edge_not(const pmg::Vert* vert) const;
    pmg::Edge* find_edge(const pmg::Vert* vert0, const pmg::Vert* vert1) const;
    pmg::Edge* shortest_edge() const;
    pmg::Edge* longest_edge() const;

    bool contains(const pmg::Edge* edge) const;
    bool contains(const pmg::Vert* vert) const;

    Face(const pmg::Edge* edge0, const pmg::Edge* edge1, const pmg::Edge* edge2);
    Face(const pmg::Vert* vert0, const pmg::Vert* vert1, const pmg::Vert* vert2);
};

} // namespace pmg
