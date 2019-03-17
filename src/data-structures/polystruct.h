// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <stddef.h>
#include <vector>


namespace polygen {

struct PolyStruct
{
    struct NodePos
    {
        double x;
        double y;
        double z;

        double& operator[](unsigned i);
        const double& operator[](unsigned i) const;
    };

    struct Facet
    {
        size_t node0;
        size_t node1;
        size_t node2;

        size_t& operator[](unsigned i);
        const size_t& operator[](unsigned i) const;
    };

    typedef size_t FacetIndex;
    typedef std::vector<FacetIndex> Crys;

    std::vector<NodePos>  nodes;
    std::vector<Facet>    facets;
    std::vector<Crys>     cryses;

    void clear();

    PolyStruct& operator=(PolyStruct&& other) noexcept;
    PolyStruct& operator=(const PolyStruct& other);

    PolyStruct(PolyStruct&& other) noexcept;
    PolyStruct(const PolyStruct& other);
    PolyStruct();
};

} // namespace polygen
