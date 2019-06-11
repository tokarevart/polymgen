#pragma once
#include <cstddef>
#include <array>
#include <vector>
#include "vec.h"


namespace spt {

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
class polytope;


template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using peak = polytope<N, Dim, Real>;

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using ridge = polytope<N, Dim, Real>;

template <std::size_t N, std::size_t Dim = 3, typename Real = typename spt::vec<Dim>::real_type>
using facet = polytope<N, Dim, Real>;


template <std::size_t N, std::size_t Dim, typename Real>
class polytope
{
public:
    static constexpr auto n = N;
    static constexpr auto dim = Dim;
    using real_type = Real;

    const std::vector<facet<N - 1, Dim, Real>*> facets;

    template <typename SubPolytope>
    bool contains(const SubPolytope* subpt) const
    {
        if constexpr (std::is_same<facet<N - 1, Dim, Real>, SubPolytope>())
            return std::find(facets.begin(), facets.end(), subpt) != facets.end();
        else
        {
            for (const auto& facet : facets)
                if (facet->contains(subpt))
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
