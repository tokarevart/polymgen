// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <array>
#include <vector>
#include "../real-type.h"


namespace psg {

struct PolyShell {
    using coordinate_t = real_t;
    using VertIdx = std::size_t;
    using FaceIdx = std::size_t;
    using Vert  = std::array<coordinate_t, 3>;
    using Face  = std::array<VertIdx, 3>;
    using Polyh = std::vector<FaceIdx>;

    std::vector<Vert>  verts;
    std::vector<Face>  faces;
    std::vector<Polyh> polyhs;

    bool empty() const;
    void clear();

    PolyShell& operator=(PolyShell&& other) noexcept;
    PolyShell& operator=(const PolyShell& other);

    PolyShell(PolyShell&& other) noexcept;
    PolyShell(const PolyShell& other);
    PolyShell();
};

} // namespace polygen
