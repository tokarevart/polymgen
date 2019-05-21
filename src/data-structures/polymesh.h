// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <array>
#include <vector>
#include "real-type.h"


namespace pmg {

struct PolyMesh
{
    using coordinate_t = real_t;
    using VertIdx = std::size_t;
    using TetrIdx = std::size_t;
    using Vert  = std::array<coordinate_t, 3>;
    using Tetr  = std::array<VertIdx, 4>;
    using Polyh = std::vector<TetrIdx>;

    std::vector<Vert>  verts;
    std::vector<Tetr>  tetrs;
    std::vector<Polyh> polyhs;

    bool empty() const;
    void clear();

    PolyMesh& operator=(PolyMesh&& other) noexcept;
    PolyMesh& operator=(const PolyMesh& other);

    PolyMesh(PolyMesh&& other) noexcept;
    PolyMesh(const PolyMesh& other);
    PolyMesh();
};

} // namespace pmg
