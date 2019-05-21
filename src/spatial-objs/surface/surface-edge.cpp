// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/surface/surface-edge.h"
#include <cmath>
#include <algorithm>
#include "helpers/spatial-algs/vec.h"

using namespace pmg;


const std::vector<Edge*>& surface::Edge::innerEdges() const
{
    return m_innerEdges;
}


const std::vector<Vert*>& surface::Edge::innerVerts() const
{
    return m_innerVerts;
}


void surface::Edge::segmentize(real_t preferredLen)
{
    std::size_t n_inner_verts = static_cast<std::size_t>(std::round(magnitude() / preferredLen)) - 1;
    if (n_inner_verts == 0)
    {
        m_innerEdges.push_back(new pmg::Edge(verts[0]->attachedVert, verts[1]->attachedVert));
        return;
    }

    vec3 dir = (*verts[1] - *verts[0]) / (n_inner_verts + 1);

    vec3 cur_pos = verts[0]->pos();
    for (std::size_t i = 0; i < n_inner_verts; i++)
    {
        cur_pos += dir;
        m_innerVerts.push_back(new pmg::Vert(cur_pos));
    }

    m_innerEdges.push_back(new pmg::Edge(verts[0]->attachedVert, m_innerVerts.front()));

    for (std::size_t i = 0; i < m_innerVerts.size() - 1; i++)
        m_innerEdges.push_back(new pmg::Edge(m_innerVerts[i], m_innerVerts[i + 1]));

    m_innerEdges.push_back(new pmg::Edge(verts[1]->attachedVert, m_innerVerts.back()));
}


real_t surface::Edge::magnitude() const
{
    return std::sqrt(sqrMagnitude());
}


real_t surface::Edge::sqrMagnitude() const
{
    vec3 buf = *verts[1] - *verts[0];
    return vec3::dot(buf, buf);
}


// TODO: replace its usage with pmg::relations content
bool surface::Edge::contains(const surface::Vert* sVert) const
{
    if (verts[0] == sVert ||
        verts[1] == sVert)
        return true;

    return false;
}


// TODO: replace its usage with pmg::relations content
bool surface::Edge::contains(const pmg::Edge* edge) const
{
    if (std::find(m_innerEdges.begin(), m_innerEdges.end(), edge) != m_innerEdges.end())
        return true;

    return false;
}


// TODO: replace its usage with pmg::relations content
bool surface::Edge::contains(const pmg::Vert* vert) const
{
    if (verts[0]->attachedVert == vert ||
        verts[1]->attachedVert == vert ||
        std::find(m_innerVerts.begin(), m_innerVerts.end(), vert) != m_innerVerts.end())
        return true;

    return false;
}




surface::Edge::Edge(const surface::Vert* vert0, const surface::Vert* vert1)
{
    verts[0] = const_cast<surface::Vert*>(vert0);
    verts[1] = const_cast<surface::Vert*>(vert1);
}
