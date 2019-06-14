#pragma once
#include <cstddef>
#include <stdexcept>
#include <array>
#include <vector>
#include "vec.h"


namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
struct polytope;

template <typename Polytope>
using aggregate = std::vector<Polytope>;


template <std::size_t N, std::size_t Dim, typename Real>
struct polytope
{
    static constexpr auto n = N;
    static constexpr auto dim = Dim;
    using real_type = Real;
    using facet_type = polytope<N - 1, Dim, Real>;

    std::vector<facet_type*> facets;

    template <typename SubPolytope>
    std::vector<SubPolytope*> all_of() const;

    bool empty() const
    {
        return facets.empty();
    }

    template <typename SubPolytope>
    bool contains(const SubPolytope* subpt) const
    {
        if constexpr (SubPolytope::n >= n)
            return false;
        else if constexpr (std::is_same<polytope<N - 1, Dim, Real>, SubPolytope>())
            return std::find(facets.begin(), facets.end(), subpt) != facets.end();
        else
        {
            for (const auto& facet : facets)
                if (facet->contains(subpt))
                    return true;

            return false;
        }
    }

    polytope(const polytope& other)
    {
        facets = other.facets;
    }
    polytope(polytope&& other) noexcept
    {
        facets = std::move(other.facets);
    }
    polytope(const std::vector<facet_type*>& facets)
        : facets(facets) {}
    polytope(std::vector<facet_type*>&& facets) noexcept
        : facets(std::move(facets)) {}
    template <std::size_t NFacets>
    polytope(const std::array<facet_type*, NFacets>& facets)
        : facets(facets.begin(), facets.end()) {}
    template <typename... Facets>
    polytope(Facets... facets)
        : facets({ facets... }) {}
};

} // namespace spt
