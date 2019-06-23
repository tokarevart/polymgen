#pragma once
#include "../polyspt/simplex.h"


namespace pmg {

template <std::size_t Dim, typename ValueType>
ValueType quality(const spt::simplex_v<3, Dim, ValueType>* simp) {
    std::array<ValueType, 4> sqr_prods = {
        static_cast<ValueType>(1.0),
        static_cast<ValueType>(1.0),
        static_cast<ValueType>(1.0),
        static_cast<ValueType>(1.0)
    };
    for (std::size_t i = 0; i < 4; i++)
        for (std::size_t j = 0; j < 4; j++)
            if (j != i)
                sqr_prods[i] *= (simp->vertices[j]->pos - simp->vertices[i]->pos).sqr_magnitude();

    ValueType max_sqr_prod = std::max({ sqr_prods[0], sqr_prods[1], sqr_prods[2], sqr_prods[3] });

    constexpr auto CONST_12_SQRT2 = static_cast<ValueType>(8.4852813742385702928101323452582);

    return CONST_12_SQRT2 * simp->volume() / std::sqrt(max_sqr_prod);
}

} // namespace pmg
