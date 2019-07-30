#pragma once
#include "mesher-base.h"


namespace pmg {

template <typename Polytope, template <typename... Args> typename Pointer>
class mesher<Polytope, spt::simplex, Pointer, spt::amount::single> {
    using vertex_type = spt::polytope<0, Polytope::dim, typename Polytope::real_type>;
    using edge_type = spt::polytope<1, Polytope::dim, typename Polytope::real_type>;
    using facet_type = spt::polytope<Polytope::n - 1, Polytope::dim, typename Polytope::real_type>;

public:
    using polytope_type = Polytope;
    using shell_type = spt::mesh_base<Pointer, facet_type>;
    using shell_mesh_unit_type = spt::simplex<Polytope::n - 1, Polytope::dim, typename Polytope::real_type>;
    using shell_mesh_type = spt::mesh_base<Pointer, shell_mesh_unit_type, spt::amount::aggregate>;
    using mesh_unit_type = spt::simplex<Polytope::n, Polytope::dim, typename Polytope::real_type>;
    using mesh_type = spt::mesh_base<Pointer, mesh_unit_type>;
    using real_type = typename polytope_type::real_type;

    void run(real_type preferred_length,
             const genparams<polytope_type>& gen_params = genparams<polytope_type>());

    mesher(shell_type shell, shell_mesh_type mesh) 
        : m_shell{ std::move(shell) }, m_shell_mesh{ std::move(mesh) } {}
    mesher(shell_type shell) : m_shell{ std::move(shell) } {
        spt::unique_mesh<facet_type> l_shell;
        if constexpr (std::is_same_v<std::unique_ptr, Pointer>)
            ;
        else
            ;
        mesher<facet_type, spt::simplex, spt::raw_ptr, spt::amount::aggregate> sh_mesher(l_shell);
        // mesh the shell...
        m_shell_mesh = std::move(/*mesh*/);
    }
    mesher(shell_mesh_type mesh) : m_shell_mesh{ std::move(mesh) } {}


private:
    real_type m_preferred_length = 0;
    genparams<polytope_type> m_gen_params;
    genparams<polytope_type> m_current_gen_params;

    shell_type m_shell;
    shell_mesh_type m_shell_mesh;
    mesh_type m_mesh;

    // methods...
};


template <typename Polytope, template <typename... Args> typename Pointer>
class mesher<Polytope, spt::simplex, Pointer, spt::amount::aggregate> {
    using vertex_type = spt::polytope<0, Polytope::dim, typename Polytope::real_type>;
    using edge_type = spt::polytope<1, Polytope::dim, typename Polytope::real_type>;
    using facet_type = spt::polytope<Polytope::n - 1, Polytope::dim, typename Polytope::real_type>;

public:
    using polytope_type = Polytope;
    using shell_type = spt::mesh_base<Pointer, facet_type, spt::amount::aggregate>;
    using shell_mesh_unit_type = spt::simplex<Polytope::n - 1, Polytope::dim, typename Polytope::real_type>;
    using shell_mesh_type = spt::raw_mesh<shell_mesh_unit_type, spt::amount::aggregate>;
    using mesh_unit_type = spt::simplex<Polytope::n, Polytope::dim, typename Polytope::real_type>;
    using mesh_type = spt::raw_mesh<mesh_unit_type, spt::amount::aggregate>;
    using real_type = typename polytope_type::real_type;

    void run(real_type preferred_length,
             const genparams<polytope_type>& gen_params = genparams<polytope_type>());

    mesher(shell_type shell, shell_mesh_type mesh)
        : m_shell{ std::move(shell) }, m_shell_mesh{ std::move(mesh) } {}
    mesher(shell_type shell) : m_shell{ std::move(shell) } {
        spt::unique_mesh<facet_type> l_shell;
        if constexpr (std::is_same_v<std::unique_ptr, Pointer>)
            ;
        else
            ;
        mesher<facet_type, spt::simplex, spt::raw_ptr, spt::amount::aggregate> sh_mesher(l_shell);
        // mesh the shell...
        m_shell_mesh = std::move(/*mesh*/);
    }
    mesher(shell_mesh_type mesh) : m_shell_mesh{ std::move(mesh) } {}


private:
    real_type m_preferred_length;
    genparams<polytope_type> m_gen_params;
    genparams<polytope_type> m_current_gen_params;

    shell_type m_shell;
    shell_mesh_type m_shell_mesh;
    mesh_type m_meshes;

    // methods...
};


template <
    template <std::size_t N, std::size_t Dim, typename Real> typename ElemType,
    std::size_t N_Minus_1, std::size_t Dim, typename Real,
    template <typename... Args> typename Pointer,
    spt::amount Amount>
mesher(const spt::mesh_base<Pointer, spt::polytope<N_Minus_1, Dim, Real>, Amount>& shell,
       const spt::mesh_base<Pointer, ElemType<N_Minus_1, Dim, Real>, spt::amount::aggregate>& shell_mesh)
    ->mesher<spt::polytope<N_Minus_1 + 1, Dim, Real>, ElemType, Pointer, Amount>;

template <
    template <std::size_t N, std::size_t Dim, typename Real> typename ElemType,
    std::size_t N_Minus_1, std::size_t Dim, typename Real,
    template <typename... Args> typename Pointer,
    spt::amount Amount>
mesher(spt::mesh_base<Pointer, spt::polytope<N_Minus_1, Dim, Real>, Amount>&& shell,
       spt::mesh_base<Pointer, ElemType<N_Minus_1, Dim, Real>, spt::amount::aggregate>&& shell_mesh)
    ->mesher<spt::polytope<N_Minus_1 + 1, Dim, Real>, ElemType, Pointer, Amount>;

} // namespace pmg