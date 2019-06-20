#pragma once
#include <cstddef>
#include <vector>
#include <array>
#include <algorithm>
#include "vertex.h"

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
