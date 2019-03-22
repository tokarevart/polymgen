// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <stddef.h>
#include <vector>
#include "real-type.h"


namespace psg {

struct PolyShell
{
    struct VertPos
    {
        real_t x;
        real_t y;
        real_t z;

        real_t& operator[](unsigned i);
        real_t  operator[](unsigned i) const;
    };

    struct Face
    {
        size_t vert0;
        size_t vert1;
        size_t vert2;

        size_t& operator[](unsigned i);
        size_t  operator[](unsigned i) const;
    };

    typedef size_t FaceIdx;
    typedef std::vector<FaceIdx> Polyhedron;

    std::vector<VertPos>    verts;
    std::vector<Face>       faces;
    std::vector<Polyhedron> polyhedrons;

    void clear();

    PolyShell& operator=(PolyShell&& other) noexcept;
    PolyShell& operator=(const PolyShell& other);

    PolyShell(PolyShell&& other) noexcept;
    PolyShell(const PolyShell& other);
    PolyShell();
};

} // namespace polygen
