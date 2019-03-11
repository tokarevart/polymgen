#include "spatial-objs/tetr.h"
#include <cmath>
#include "helpers/spatial-algs/vec.h"

using pmg::Tetr;




double Tetr::computeVolume() const
{
    tva::Vec v0 = verts[0]->pos();
    return 0.16666666666666666 * std::abs(tva::Vec::mixed(
        verts[1]->pos() - v0,
        verts[2]->pos() - v0,
        verts[3]->pos() - v0));
}

double Tetr::computeQuality() const
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


Tetr::Tetr(const Vertex* vert0, const Vertex* vert1, const Vertex* vert2, const Vertex* vert3)
{
    verts[0] = const_cast<Vertex*>(vert0);
    verts[1] = const_cast<Vertex*>(vert1);
    verts[2] = const_cast<Vertex*>(vert2);
    verts[3] = const_cast<Vertex*>(vert3);
}
