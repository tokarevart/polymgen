#pragma once
#include <cstddef>
#include <vector>
#include <array>
#include <algorithm>
#include "edge.h"

namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
struct simplex
{
    static constexpr auto n = N;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using facet_type = simplex<N - 1, Dim, Real>;

    std::array<facet_type*, N + 1> facets;
   
    template <typename SubPolytope>
    std::vector<SubPolytope*> all_of() const;

    bool empty() const
    {
        return std::none_of(facets.begin(), facets.end(), [](auto x) -> bool { return x; });
    }

    template <typename SubSimplex>
    bool contains(const SubSimplex* subsimp) const
    {
        if constexpr (SubSimplex::n >= n)
            return false;
        else if constexpr (std::is_same<facet_type, SubSimplex>())
            return std::find(facets.begin(), facets.end(), subsimp) != facets.end();
        else
        {
            for (const auto& facet : facets)
                if (facet->contains(subsimp))
                    return true;

            return false;
        }
    }

    simplex(const simplex& other)
    {
        facets = other.facets;
    }
    simplex(const std::array<facet_type*, N + 1>& facets)
        : facets(facets) {}
    template <typename... Facets>
    simplex(Facets... facets)
        : facets({ facets... }) {}
};

template <std::size_t Dim, typename Real>
struct simplex<2, Dim, Real>
{
    static constexpr auto n = 2;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using facet_type = spt::edge<Dim, Real>;

    std::array<facet_type*, 3> facets = { nullptr, };

    bool empty() const
    {
        return !facets[0] && !facets[1] && !facets[2];
    }

    template <typename SubPolytope>
    auto all_of() const
    {
        static_assert(SubPolytope::n == 1 || SubPolytope::n == 0);
        if constexpr (SubPolytope::n == 1)
            return facets;
        else
        {
            std::array<SubPolytope*, 4> res = { nullptr, };
            //
        }
    }

    bool contains(const facet_type* subpt) const
    {
        return facets[0] == subpt || facets[1] == subpt || facets[2] == subpt;
    }

    simplex(const simplex& poly)
    {
        facets = poly.facets;
    }
    simplex(const std::array<facet_type*, 2>& facets)
        : facets(facets) {}
    simplex(const facet_type* f0, const facet_type* f1, const facet_type* f2)
        : facets({ f0, f1, f2 }) {}
};

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
struct simplex_v
{
    static constexpr auto n = N;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using vertex_type = spt::vertex<Dim, Real>;

    std::array<vertex_type*, N + 1> vertices;

    bool empty() const
    {
        return std::none_of(vertices.begin(), vertices.end(), [](auto x) -> bool { return x; });
    }

    bool contains(const vertex_type* vert) const
    {
        return std::find(vertices.begin(), vertices.end(), vert) != vertices.end();
    }

    simplex_v(const simplex_v& other)
    {
        vertices = other.vertices;
    }
    simplex_v(const std::array<vertex_type*, N + 1>& verts)
        : vertices(vertices) {}
    template <typename... Vertices>
    simplex_v(Vertices... verts)
        : vertices({ verts... }) {}
};

} // namespace spt
