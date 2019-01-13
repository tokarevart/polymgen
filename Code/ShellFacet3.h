#pragma once
#include <vector>
#include <memory>
#include "Definitions.h"
#include "Inclusions.h"

class ShellFacet3
{
public:
    ShellEdge3* edges[3];

    ShellVertex3* findVertexNotIncludedInEdge(const ShellEdge3& edge) const;

    const bool contains(const ShellEdge3& edge) const;
    const bool contains(const ShellVertex3& vertex) const;


    ShellFacet3();
    ShellFacet3(
        ShellEdge3& edge0, 
        ShellEdge3& edge1, 
        ShellEdge3& edge2);
    ~ShellFacet3();
};