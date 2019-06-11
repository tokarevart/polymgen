#pragma once
#include <cstddef>
#include <array>
#include "vec.h"

namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
class simplex
{
public:
    static constexpr auto n = N;
    static constexpr auto dim = Dim;
    using real_type = Real;

    const std::array<simplex<N - 1, Dim, Real>*, N + 1> facets;

    template <typename SubPolytope>
    std::array<SubPolytope*, /*How many? I feel there is a formula*/> all_of() const;

    template <typename SubSimplex>
    bool contains(const SubSimplex* subsimp) const
    {
        if constexpr (SubSimplex::n >= n)
            return false;
        else if constexpr (std::is_same<simplex<N - 1, Dim, Real>, SubSimplex>())
            return std::find(facets.begin(), facets.end(), subsimp) != facets.end();
        else
        {
            for (const auto& facet : facets)
                if (facet->contains(subsimp))
                    return true;

            return false;
        }
    }

    polytope() = delete;
    polytope(const std::vector<facet<N - 1, Dim, Real>*>& facets)
        : facets(facets) {}
    polytope(std::vector<facet<N - 1, Dim, Real>*>&& facets) noexcept
        : facets(std::move(facets)) {}
    template <std::size_t NFacets>
    polytope(const std::array<facet<N - 1, Dim, Real>*, NFacets>& facets)
        : facets(facets.begin(), facets.end()) {}
    template <typename... Facets>
    polytope(Facets... facets)
        : facets({ facets... }) {}
};

} // namespace spt
