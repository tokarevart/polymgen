#include "Simplex3.h"
#include <cmath>
#include "helpers/spatialalgs/Vec3.h"

namespace tva { struct Vec3; }


double Simplex3::computeVolume() const
{
    tva::Vec3 v0 = verts[0]->getPos();
    return 0.16666666666666666 * std::abs(tva::Vec3::mixedProduct(
        verts[1]->getPos() - v0, 
        verts[2]->getPos() - v0, 
        verts[3]->getPos() - v0));
}

double Simplex3::computeQuality() const
{
    double sqr_prods[4];
    for (int i = 0; i < 4; i++)
    {
        sqr_prods[i] = 1.0;
        for (int j = 0; j < 4; j++)
            if (j != i)
                sqr_prods[i] *= (*verts[j] - *verts[i]).sqrMagnitude();
    }
    double max_sqr_prod = std::max({ sqr_prods[0], sqr_prods[1], sqr_prods[2], sqr_prods[3] });

    return 8.48528137423857 * computeVolume() / std::sqrt(max_sqr_prod);
}

Simplex3::Simplex3() {}

Simplex3::Simplex3(const Vertex3* vert0, const Vertex3* vert1, const Vertex3* vert2, const Vertex3* vert3)
{
    verts[0] = const_cast<Vertex3*>(vert0);
    verts[1] = const_cast<Vertex3*>(vert1);
    verts[2] = const_cast<Vertex3*>(vert2);
    verts[3] = const_cast<Vertex3*>(vert3);
}

Simplex3::~Simplex3() {}
