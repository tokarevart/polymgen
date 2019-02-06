#pragma once
#include <stddef.h>
#include <vector>


struct PolyStruct
{
    struct NodePos
    {
        double x;
        double y;
        double z;

        double& operator[](unsigned i)
        {
            switch (i)
            {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: return z;
            }
        }

        const double& operator[](unsigned i) const
        {
            return const_cast<NodePos&>(*this)[i];
        }
    };

    struct Facet
    {
        size_t node0;
        size_t node1;
        size_t node2;

        size_t& operator[](unsigned i)
        {
            switch (i)
            {
            case 0: return node0;
            case 1: return node1;
            case 2: return node2;
            default: return node2;
            }
        }

        const size_t& operator[](unsigned i) const
        {
            return const_cast<Facet&>(*this)[i];
        }
    };

    typedef size_t FacetIndex;
    typedef std::vector<FacetIndex> Crys;

    std::vector<NodePos>  nodes;
    std::vector<Facet>    facets;
    std::vector<Crys>     cryses;

    void clear();

    PolyStruct& operator=(PolyStruct&& right) noexcept;
    PolyStruct& operator=(const PolyStruct& right);

    PolyStruct(PolyStruct&& polyStruct) noexcept;
    PolyStruct(const PolyStruct& polyStruct);
    PolyStruct() {}
    ~PolyStruct() {}
};
