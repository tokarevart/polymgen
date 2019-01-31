#pragma once

struct PolyStruct
{
    size_t  nNodes;
    double* nodesPositions; // { x0, y0, z0, x1, y1, z1 ... }
    size_t  nFacets;
    size_t* facets;         // { facet0_node0, facet0_node1, facet0_node2, facet1_node0, facet1_node1, facet1_node2 ... }
    size_t  nCryses;
    size_t* nCrysesFacets;  // Number of facets in each crystallite's shell.
    size_t* cryses;         // { crys0_facet0, crys0_facet1, ... crys1_facet0, crys1_facet1 ... }

    PolyStruct();
    ~PolyStruct();
};