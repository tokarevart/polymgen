#pragma once
#include "mesher-base.h"


namespace pmg {

template <typename Polytope>
class mesher<Polytope, spt::elem_shape::simplex> {
    using vertex_type = spt::polytope<0, Polytope::dim, typename Polytope::value_type>;
    using edge_type = spt::polytope<1, Polytope::dim, typename Polytope::value_type>;
    using facet_type = spt::polytope<Polytope::n - 1, Polytope::dim, typename Polytope::value_type>;

public:
    using polytope_type = Polytope;
    using shell_type = polytope_type;
    using shell_mesh_type = spt::mesh<spt::aggregate<facet_type*>, spt::elem_shape::simplex>;
    using mesh_type = spt::mesh<polytope_type, spt::elem_shape::simplex>;
    using value_type = typename polytope_type::value_type;

    void run(value_type preferred_length,
             const genparams<polytope_type>& gen_params = genparams<polytope_type>());


    mesher(const shell_type& shell, const shell_mesh_type& mesh) {
        m_shell = shell;
        m_shell_mesh = mesh;
    }
    mesher(shell_type&& shell, shell_mesh_type&& mesh) noexcept {
        m_shell = std::move(shell);
        m_shell_mesh = std::move(mesh);
    }
    mesher(const shell_type& shell) {
        m_shell = shell;
        mesher<spt::aggregate<facet_type*>, spt::elem_shape::simplex> sh_mesher(shell);
        // mesh the shell...
        m_shell_mesh = std::move(/*mesh*/);
    }
    mesher(shell_type&& shell) noexcept {
        m_shell = std::move(shell);
        mesher<spt::aggregate<facet_type*>, spt::elem_shape::simplex> sh_mesher(shell);
        // mesh the shell...
        m_shell_mesh = std::move(/*mesh*/);
    }
    mesher(const shell_mesh_type& mesh) {
        m_shell_mesh = mesh;
    }


private:
    value_type m_preferred_length;
    genparams<polytope_type> m_gen_params;
    genparams<polytope_type> m_current_gen_params;

    shell_type m_shell;
    shell_mesh_type m_shell_mesh;
    mesh_type m_mesh;

    // methods...
};


template <typename Polytope>
class mesher<spt::composition<Polytope>, spt::elem_shape::simplex> {
    using vertex_type = spt::polytope<0, Polytope::dim, typename Polytope::value_type>;
    using edge_type = spt::polytope<1, Polytope::dim, typename Polytope::value_type>;
    using facet_type = spt::polytope<N - 1, Polytope::dim, typename Polytope::value_type>;

public:
    using polytope_type = Polytope;
    using shell_type = spt::aggregate<polytope_type>;
    using shell_mesh_type = spt::mesh<spt::aggregate<facet_type>, spt::elem_shape::simplex>;
    using mesh_type = spt::mesh<spt::aggregate<polytope_type>, spt::elem_shape::simplex>;
    using value_type = typename polytope_type::value_type;

    void run(value_type preferred_length,
             const genparams<polytope_type>& gen_params = genparams<polytope_type>());


    mesher(const shell_type& shell, const shell_mesh_type& mesh) {
        m_shell = shell;
        m_shell_mesh = mesh;
    }
    mesher(shell_type&& shell, shell_mesh_type&& mesh) noexcept {
        m_shell = std::move(shell);
        m_shell_mesh = std::move(mesh);
    }
    mesher(const shell_type& shell) {
        m_shell = shell;
        mesher<spt::aggregate<facet_type>, spt::elem_shape::simplex> sh_mesher(shell);
        // mesh the shell...
        m_shell_mesh = std::move(/*mesh*/);
    }
    mesher(shell_type&& shell) noexcept {
        m_shell = std::move(shell);
        mesher<spt::aggregate<facet_type>, spt::elem_shape::simplex> sh_mesher(shell);
        // mesh the shell...
        m_shell_mesh = std::move(/*mesh*/);
    }
    mesher(const shell_mesh_type& mesh) {
        m_shell_mesh = mesh;
    }


private:
    value_type m_preferred_length;
    genparams<polytope_type> m_gen_params;
    genparams<polytope_type> m_current_gen_params;

    shell_type m_shell;
    shell_mesh_type m_shell_mesh;
    mesh_type m_meshes;

    // methods...
};

} // namespace pmg