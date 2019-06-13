#pragma once
#include "polytope.h"
#include "mesh.h"
#include "genparams.h"


namespace pmg {

template <typename Polytope>
class mesher
{
    using vertex_type = spt::polytope<0, Polytope::dim, typename Polytope::real_type>;
    using edge_type = spt::polytope<1, Polytope::dim, typename Polytope::real_type>;
    using facet_type = spt::polytope<Polytope::n - 1, Polytope::dim, typename Polytope::real_type>;
    using shell_mesh_type = pmg::mesh<spt::aggregate<facet_type>>;
    using mesh_type = pmg::mesh<facet_type>;

public:
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;

    void run(typename real_type preferred_length, 
        const genparams<polytope_type>& gen_params = genparams<polytope_type>());

    mesher(const shell_mesh_type* mesh)
    {
        copy_mesh(mesh);
    }

    mesher(const shell_mesh_type* mesh,
        typename real_type preferred_length,
        const genparams<polytope_type>& gen_params = genparams<polytope_type>())
    {
        copy_mesh(mesh);
        run(preferred_length, gen_params);
    }


private:
    real_type m_preferred_length;
    genparams<polytope_type> m_gen_params;
    genparams<polytope_type> m_current_gen_params;

    shell_mesh_type m_shell; // or 
    mesh_type m_mesh;

    void copy_mesh(const raw_mesh_type& mesh);
};


template <typename Polytope>
class mesher<spt::aggregate<Polytope>>
{
    using vertex_type = spt::polytope<0, Polytope::dim, typename Polytope::real_type>;
    using edge_type = spt::polytope<1, Polytope::dim, typename Polytope::real_type>;
    using facet_type = spt::polytope<Polytope::n - 1, Polytope::dim, typename Polytope::real_type>;
    using mesh_type = pmg::mesh<facet_type>;

public:
    using polytope_type = Polytope;
    using real_type = typename polytope_type::real_type;

    void run(typename real_type preferred_length,
        const genparams<polytope_type>& gen_params = genparams<polytope_type>());

    mesher(const mesh_type* mesh)
    {
        copy_mesh(mesh);
    }

    mesher(const mesh_type* mesh,
        typename real_type preferred_length,
        const genparams<polytope_type>& gen_params = genparams<polytope_type>())
    {
        copy_mesh(mesh);
        run(preferred_length, gen_params);
    }


private:
    real_type m_preferred_length;
    genparams<polytope_type> m_gen_params;
    genparams<polytope_type> m_current_gen_params;

    mesh_type m_shell; // or 
    mesh_type m_mesh;

    void copy_mesh(const raw_mesh_type& mesh);
};

} // namespace pmg