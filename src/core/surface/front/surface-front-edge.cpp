// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "edge.h"
#include "../../../helpers/spatial/algs.h"

using namespace pmg::surface;
using vec3 = spt::vec<3, real_t>;


vec3 front::Edge::compute_normal() {
    surface::Edge* related_sedge = m_related_surface_face->find_surface_edge_containing(x);
    vec3 opp_v_pos = m_related_surface_face->find_vert_not(related_sedge)->pos();
    vec3 p0 = related_sedge->verts[0]->pos();
    vec3 p1 = related_sedge->verts[1]->pos();

    return normal = (opp_v_pos - spt::project(opp_v_pos, p0, p1)).normalize();
}


vec3 front::Edge::center() {
    return (x->verts[0]->pos() + x->verts[1]->pos()) * static_cast<real_t>(0.5);
}




front::Edge::Edge(const pmg::surface::Face* relatedSurfaceFace, const pmg::Edge* edge)
    : x(const_cast<pmg::Edge*>(edge)), m_related_surface_face(const_cast<surface::Face*>(relatedSurfaceFace)) {}
