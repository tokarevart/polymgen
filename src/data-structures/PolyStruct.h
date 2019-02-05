#pragma once
#include <stddef.h>
#include <vector>


//struct PolyStruct
//{
//    size_t  nNodes;
//    double* nodesPositions; // { x0, y0, z0, x1, y1, z1 ... }
//    size_t  nFacets;
//    size_t* facets;         // { facet0_node0, facet0_node1, facet0_node2, facet1_node0, facet1_node1, facet1_node2 ... }
//    size_t  nCryses;
//    size_t* nCrysesFacets;  // Number of facets in each crystallite's shell.
//    size_t* cryses;         // { crys0_facet0, crys0_facet1, ... crys1_facet0, crys1_facet1 ... }
//};

struct PolyStruct
{
    struct NodePos
    {
        double x;
        double y;
        double z;

        double& operator[](int i)
        {
            switch (i)
            {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: return z;
            }
        }

        const double& operator()(int i) const
        {
            switch (i)
            {
            case 0: return x;
            case 1: return y;
            case 2: return z;
            default: return z;
            }
        }
    };

    struct Facet
    {
        size_t node0;
        size_t node1;
        size_t node2;

        size_t& operator[](int i)
        {
            switch (i)
            {
            case 0: return node0;
            case 1: return node1;
            case 2: return node2;
            default: return node2;
            }
        }

        const size_t& operator()(int i) const
        {
            switch (i)
            {
            case 0: return node0;
            case 1: return node1;
            case 2: return node2;
            default: return node2;
            }
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
