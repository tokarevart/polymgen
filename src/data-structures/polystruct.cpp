// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "data-structures/polystruct.h"

using namespace polygen;




real_t& PolyStruct::NodePos::operator[](unsigned i)
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return z;
    }
}


real_t PolyStruct::NodePos::operator[](unsigned i) const
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return z;
    }
}




size_t& PolyStruct::Facet::operator[](unsigned i)
{
    switch (i)
    {
    case 0: return node0;
    case 1: return node1;
    case 2: return node2;
    default: return node2;
    }
}


size_t PolyStruct::Facet::operator[](unsigned i) const
{
    switch (i)
    {
    case 0: return node0;
    case 1: return node1;
    case 2: return node2;
    default: return node2;
    }
}




void PolyStruct::PolyStruct::clear()
{
    nodes.clear();
    facets.clear();
    cryses.clear();
}




PolyStruct& PolyStruct::PolyStruct::operator=(PolyStruct&& other) noexcept
{
    nodes  = std::move(other.nodes);
    facets = std::move(other.facets);
    cryses = std::move(other.cryses);
    return *this;
}


PolyStruct& PolyStruct::PolyStruct::operator=(const PolyStruct& other)
{
    nodes  = other.nodes;
    facets = other.facets;
    cryses = other.cryses;
    return *this;
}




PolyStruct::PolyStruct::PolyStruct(PolyStruct&& other) noexcept
    : nodes(std::move(other.nodes)), facets(std::move(other.facets)), cryses(std::move(other.cryses)) {}


PolyStruct::PolyStruct::PolyStruct(const PolyStruct& other)
{
    nodes  = other.nodes;
    facets = other.facets;
    cryses = other.cryses;
}


PolyStruct::PolyStruct::PolyStruct() {}
