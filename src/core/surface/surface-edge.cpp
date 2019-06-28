// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "edge.h"
#include <cmath>
#include <algorithm>
#include "../../helpers/spatial/algs.h"

using namespace pmg;
using vec3 = spt::vec<3, real_t>;


const std::vector<Edge*>& surface::Edge::inner_edges() const {
    return m_inner_edges;
}


const std::vector<Vert*>& surface::Edge::inner_verts() const {
    return m_inner_verts;
}


void surface::Edge::segmentize(real_t preferredLen) {
    std::size_t n_inner_verts = static_cast<std::size_t>(std::round(magnitude() / preferredLen)) - 1;
    if (n_inner_verts == 0) {
        m_inner_edges.push_back(new pmg::Edge(verts[0]->attached_vert, verts[1]->attached_vert));
        return;
    }

    vec3 dir = (verts[1]->pos() - verts[0]->pos()) / (n_inner_verts + 1);

    vec3 cur_pos = verts[0]->pos();
    for (std::size_t i = 0; i < n_inner_verts; i++) {
        cur_pos += dir;
        m_inner_verts.push_back(new pmg::Vert(cur_pos));
    }

    m_inner_edges.push_back(new pmg::Edge(verts[0]->attached_vert, m_inner_verts.front()));

    for (std::size_t i = 0; i < m_inner_verts.size() - 1; i++)
        m_inner_edges.push_back(new pmg::Edge(m_inner_verts[i], m_inner_verts[i + 1]));

    m_inner_edges.push_back(new pmg::Edge(verts[1]->attached_vert, m_inner_verts.back()));
}


real_t surface::Edge::magnitude() const {
    return std::sqrt(sqr_magnitude());
}


real_t surface::Edge::sqr_magnitude() const {
    vec3 buf = verts[1]->pos() - verts[0]->pos();
    return spt::dot(buf, buf);
}


bool surface::Edge::contains(const surface::Vert* sVert) const {
    if (verts[0] == sVert ||
        verts[1] == sVert)
        return true;

    return false;
}


bool surface::Edge::contains(const pmg::Edge* edge) const {
    if (std::find(m_inner_edges.begin(), m_inner_edges.end(), edge) != m_inner_edges.end())
        return true;

    return false;
}


bool surface::Edge::contains(const pmg::Vert* vert) const {
    if (verts[0]->attached_vert == vert ||
        verts[1]->attached_vert == vert ||
        std::find(m_inner_verts.begin(), m_inner_verts.end(), vert) != m_inner_verts.end())
        return true;

    return false;
}




surface::Edge::Edge(const surface::Vert* vert0, const surface::Vert* vert1) {
    verts[0] = const_cast<surface::Vert*>(vert0);
    verts[1] = const_cast<surface::Vert*>(vert1);
}
