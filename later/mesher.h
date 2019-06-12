#pragma once
#include "polytope.h"
#include "mesh.h"
#include "genparams.h"


namespace pmg {

template <typename Polytope>
class mesher
{
public:
    using polytope_type = Polytope;

    mesher(
        const pmg::raw_mesh<spt::polytope<
            polytope_type::n - 1,
            polytope_type::dim,
            polytope_type::real_type>>& raw_mesh,
        polytope_type::real_type preferred_length,
        genparams<Polytope> gen_params = genparams<Polytope>())
        : m_preferred_length(preferred_length), m_gen_params(gen_params)
    {
        run();
    }


private:
    const polytope_type::real_type m_preferred_length;
    const genparams<Polytope> m_gen_params;
    genparams<Polytope> m_current_gen_params;

    void run(polytope_type::real_type preferred_length, genparams<Polytope> gen_params = genparams<Polytope>());
};

template <typename Polytope>
class mesher<std::vector<Polytope*>>
{
public:
    using polytope_type = Polytope;

    mesher(
        const pmg::raw_mesh<std::vector<spt::polytope<
            polytope_type::n - 1, 
            polytope_type::dim, 
            polytope_type::real_type>>>& raw_mesh,
        polytope_type::real_type preferred_length,
        genparams<Polytope> gen_params = genparams<Polytope>())
        : m_preferred_length(preferred_length), m_gen_params(gen_params) 
    {
        run();
    }


private:
    const polytope_type::real_type m_preferred_length;
    const genparams<Polytope> m_gen_params;
    genparams<Polytope> m_current_gen_params;

    void run(polytope_type::real_type preferred_length, genparams<Polytope> gen_params = genparams<Polytope>());
};

} // namespace pmg