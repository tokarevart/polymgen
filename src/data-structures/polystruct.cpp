// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#include "data-structures/polystruct.h"

using namespace polygen;
typedef PolyStruct::NodePos NodePos;
typedef PolyStruct::Facet Facet;
typedef PolyStruct::Crys Crys;
typedef PolyStruct::FacetIndex FacetIndex;




double& NodePos::operator[](unsigned i)
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return z;
    }
}

const double& NodePos::operator[](unsigned i) const
{
    return const_cast<NodePos&>(*this)[i];
}




size_t& Facet::operator[](unsigned i)
{
    switch (i)
    {
    case 0: return node0;
    case 1: return node1;
    case 2: return node2;
    default: return node2;
    }
}

const size_t& Facet::operator[](unsigned i) const
{
    return const_cast<Facet&>(*this)[i];
}




void PolyStruct::clear()
{
    nodes.clear();
    facets.clear();
    cryses.clear();
}




PolyStruct& PolyStruct::operator=(PolyStruct&& other) noexcept
{
    nodes  = std::move(other.nodes);
    facets = std::move(other.facets);
    cryses = std::move(other.cryses);
    return *this;
}


PolyStruct& PolyStruct::operator=(const PolyStruct& other)
{
    nodes  = other.nodes;
    facets = other.facets;
    cryses = other.cryses;
    return *this;
}




PolyStruct::PolyStruct(PolyStruct&& other) noexcept
    : nodes(std::move(other.nodes)), facets(std::move(other.facets)), cryses(std::move(other.cryses))
{
    int a = 24;
    a = a + 23;
}

PolyStruct::PolyStruct(const PolyStruct& other)
{
    nodes  = other.nodes;
    facets = other.facets;
    cryses = other.cryses;
}

PolyStruct::PolyStruct() {}
