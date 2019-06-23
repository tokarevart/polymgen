#pragma once


namespace pmg {

template <typename Polytope>
struct genparams {
    using polytope_type = Polytope;
    using facet_type = typename polytope_type::facet_type;

    genparams<facet_type> facet_genparams = genparams<facet_type>();

    // mesh generation parameters...
};

} // namespace pmg
