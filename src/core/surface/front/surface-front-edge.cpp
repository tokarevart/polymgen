// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "core/surface/front/edge.h"
#include "helpers/spatial/algs.h"

using namespace pmg::surface;
using spt::vec3;


vec3 front::Edge::computeNormal()
{
    surface::Edge* related_sedge = m_relatedSurfaceFace->findSurfaceEdgeContaining(edge);
    vec3 opp_v_pos = m_relatedSurfaceFace->findVertNot(related_sedge)->pos();
    vec3 p0 = related_sedge->verts[0]->pos();
    vec3 p1 = related_sedge->verts[1]->pos();

    return normal = (opp_v_pos - spt::algs::project(opp_v_pos, p0, p1)).normalize();
}


vec3 front::Edge::computeCenter()
{
    return (edge->verts[0]->pos() + edge->verts[1]->pos()) * static_cast<real_t>(0.5);
}




front::Edge::Edge(const surface::Face* relatedSurfaceFace, const pmg::Edge* edge)
    : edge(const_cast<pmg::Edge*>(edge)), m_relatedSurfaceFace(const_cast<surface::Face*>(relatedSurfaceFace)) {}
