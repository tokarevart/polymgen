#pragma once
#include "ShellEdge3.h"
#include "ShellVertex3.h"

class ShellEdge3;
class ShellVertex3;


class ShellFacet3
{
public:
    ShellEdge3* edges[3];

    ShellVertex3* findVertNotIncludedInEdge(const ShellEdge3* edge) const;

    bool contains(const ShellEdge3* edge) const;
    bool contains(const ShellVertex3* vert) const;

    ShellFacet3();
    ShellFacet3(const ShellEdge3* edge0, const ShellEdge3* edge1, const ShellEdge3* edge2);
    ~ShellFacet3();
};