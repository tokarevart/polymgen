#include "ShellFacet3.h"


ShellVertex3* ShellFacet3::findVertNotIncludedInEdge(const ShellEdge3* edge) const
{
    for (auto& facet_edge : edges)
    {
        if (edge != facet_edge)
        {
            if (!edge->contains(facet_edge->verts[0]))
                return facet_edge->verts[0];
            else
                return facet_edge->verts[1];
        }
    }

    return nullptr;
}

bool ShellFacet3::contains(const ShellEdge3* edge) const
{
    for (auto& edge_ : edges)
        if (edge_ == edge)
            return true;

    return false;
}

bool ShellFacet3::contains(const ShellVertex3* vert) const
{
    for (auto& edge : edges)
        if (edge->contains(vert))
            return true;

    return false;
}

ShellFacet3::ShellFacet3() {}

ShellFacet3::ShellFacet3(
    const ShellEdge3* edge0,
    const ShellEdge3* edge1,
    const ShellEdge3* edge2)
{
    edges[0] = const_cast<ShellEdge3*>(edge0);
    edges[1] = const_cast<ShellEdge3*>(edge1);
    edges[2] = const_cast<ShellEdge3*>(edge2);
}

ShellFacet3::~ShellFacet3() {}
