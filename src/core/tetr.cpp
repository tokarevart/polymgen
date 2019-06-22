// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "tetr.h"
#include <cmath>
#include "../helpers/spatial/vec.h"

using spt::vec3;


real_t pmg::Tetr::computeVolume() const
{
    vec3 v0 = verts[0]->pos();
    return static_cast<real_t>(1.0 / 6.0) * std::abs(vec3::mixed(
        verts[1]->pos() - v0,
        verts[2]->pos() - v0,
        verts[3]->pos() - v0));
}

real_t pmg::Tetr::computeQuality() const
{
    std::array<real_t, 4> sqr_prods = 
    {
        static_cast<real_t>(1.0),
        static_cast<real_t>(1.0),
        static_cast<real_t>(1.0),
        static_cast<real_t>(1.0)
    };
    for (std::size_t i = 0; i < 4; i++)
        for (std::size_t j = 0; j < 4; j++)
            if (j != i)
                sqr_prods[i] *= (verts[j]->pos() - verts[i]->pos()).sqrMagnitude();

    real_t max_sqr_prod = std::max({ sqr_prods[0], sqr_prods[1], sqr_prods[2], sqr_prods[3] });

    constexpr auto CONST_12_SQRT2 = static_cast<real_t>(8.4852813742385702928101323452582);

    return CONST_12_SQRT2 * computeVolume() / std::sqrt(max_sqr_prod);
}




pmg::Tetr::Tetr(const Vert* vert0, const Vert* vert1, const Vert* vert2, const Vert* vert3)
{
    verts[0] = const_cast<Vert*>(vert0);
    verts[1] = const_cast<Vert*>(vert1);
    verts[2] = const_cast<Vert*>(vert2);
    verts[3] = const_cast<Vert*>(vert3);
}
