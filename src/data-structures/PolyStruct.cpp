#include "PolyStruct.h"




void PolyStruct::clear()
{
    nodes.clear();
    facets.clear();
    cryses.clear();
}




PolyStruct& PolyStruct::operator=(PolyStruct&& right) noexcept
{
    nodes  = std::move(right.nodes);
    facets = std::move(right.facets);
    cryses = std::move(right.cryses);
    return *this;
}


PolyStruct& PolyStruct::operator=(const PolyStruct& right)
{
    nodes  = right.nodes;
    facets = right.facets;
    cryses = right.cryses;
    return *this;
}




PolyStruct::PolyStruct(PolyStruct&& polyStruct) noexcept
    : nodes(std::move(polyStruct.nodes)), facets(std::move(polyStruct.facets)), cryses(std::move(polyStruct.cryses))
{
    int a = 24;
    a = a + 23;
}

PolyStruct::PolyStruct(const PolyStruct& polyStruct)
{
    nodes  = polyStruct.nodes;
    facets = polyStruct.facets;
    cryses = polyStruct.cryses;
}
