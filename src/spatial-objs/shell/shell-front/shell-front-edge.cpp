// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/shell/shell-front/shell-front-edge.h"
#include "helpers/spatial-algs/spatial-algs.h"

using FrPlEdge = pmg::shell::front::Edge;




Vec FrPlEdge::computeNormal()
{
    shell::Edge* related_sedge = m_relatedShellFace->findShellEdgeContaining(edge);
    Vec opp_v_pos = m_relatedShellFace->findVertNot(related_sedge)->pos();
    Vec p0 = related_sedge->verts[0]->pos();
    Vec p1 = related_sedge->verts[1]->pos();

    return normal = (opp_v_pos - spatalgs::project(opp_v_pos, p0, p1)).normalize();
}


Vec FrPlEdge::computeCenter()
{
    return static_cast<real_t>(0.5) * (edge->verts[0]->pos() + edge->verts[1]->pos());
}




bool FrPlEdge::isAdj(const FrPlEdge* edge0, const FrPlEdge* edge1)
{
    if (edge0->edge->verts[0] == edge1->edge->verts[0] ||
        edge0->edge->verts[0] == edge1->edge->verts[1] ||
        edge0->edge->verts[1] == edge1->edge->verts[0] ||
        edge0->edge->verts[1] == edge1->edge->verts[1])
        return true;

    return false;
}




FrPlEdge::Edge(const shell::Face* relatedShellFace, const pmg::Edge* edge)
    : edge(const_cast<pmg::Edge*>(edge)), m_relatedShellFace(const_cast<shell::Face*>(relatedShellFace)) {}
