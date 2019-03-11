#pragma once
#include <algorithm>
#include "Vertex3.h"

class Vertex3;


class Simplex3
{
public:
    Vertex3* verts[4];

    double computeVolume() const;
    double computeQuality() const;

    Simplex3();
    Simplex3(
        const Vertex3* vert0,
        const Vertex3* vert1,
        const Vertex3* vert2,
        const Vertex3* vert3);
    ~Simplex3();
};