// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/shell/shell-front/shell-front-edge.h"
#include "helpers/spatial-algs/spatial-algs.h"

using namespace pmg::shell;


vec3 front::Edge::computeNormal()
{
    shell::Edge* related_sedge = m_relatedShellFace->findShellEdgeContaining(edge);
    vec3 opp_v_pos = m_relatedShellFace->findVertNot(related_sedge)->pos();
    vec3 p0 = related_sedge->verts[0]->pos();
    vec3 p1 = related_sedge->verts[1]->pos();

    return normal = (opp_v_pos - spatalgs::project(opp_v_pos, p0, p1)).normalize();
}


vec3 front::Edge::computeCenter()
{
    return static_cast<real_t>(0.5) * (edge->verts[0]->pos() + edge->verts[1]->pos());
}




bool front::Edge::isAdj(const front::Edge* edge0, const front::Edge* edge1)
{
    if (edge0->edge->verts[0] == edge1->edge->verts[0] ||
        edge0->edge->verts[0] == edge1->edge->verts[1] ||
        edge0->edge->verts[1] == edge1->edge->verts[0] ||
        edge0->edge->verts[1] == edge1->edge->verts[1])
        return true;

    return false;
}




front::Edge::Edge(const shell::Face* relatedShellFace, const pmg::Edge* edge)
    : edge(const_cast<pmg::Edge*>(edge)), m_relatedShellFace(const_cast<shell::Face*>(relatedShellFace)) {}
