// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/shell/shell-edge.h"
#include <cmath>
#include <algorithm>
#include "helpers/spatial-algs/vec.h"

using namespace pmg;




const std::vector<Edge*>&shell::Edge::innerEdges() const
{
    return m_innerEdges;
}


const std::vector<Vert*>&shell::Edge::innerVerts() const
{
    return m_innerVerts;
}


void shell::Edge::segmentize(real_t preferredLen)
{
    size_t n_inner_verts = static_cast<size_t>(std::round(magnitude() / preferredLen)) - 1;
    if (n_inner_verts == 0)
    {
        m_innerEdges.push_back(new pmg::Edge(verts[0]->attachedVert, verts[1]->attachedVert));
        return;
    }

    vec3 dir = (*verts[1] - *verts[0]) / (n_inner_verts + 1);

    vec3 cur_pos = verts[0]->pos();
    for (size_t i = 0; i < n_inner_verts; i++)
    {
        cur_pos += dir;
        m_innerVerts.push_back(new pmg::Vert(cur_pos));
    }

    m_innerEdges.push_back(new pmg::Edge(verts[0]->attachedVert, m_innerVerts.front()));

    for (size_t i = 0; i < m_innerVerts.size() - 1; i++)
        m_innerEdges.push_back(new pmg::Edge(m_innerVerts[i], m_innerVerts[i + 1]));

    m_innerEdges.push_back(new pmg::Edge(verts[1]->attachedVert, m_innerVerts.back()));
}




real_t shell::Edge::magnitude() const
{
    return std::sqrt(sqrMagnitude());
}


real_t shell::Edge::sqrMagnitude() const
{
    vec3 buf = *verts[1] - *verts[0];
    return vec3::dot(buf, buf);
}




bool shell::Edge::contains(const shell::Vert* sVert) const
{
    if (verts[0] == sVert ||
        verts[1] == sVert)
        return true;

    return false;
}


bool shell::Edge::contains(const pmg::Edge* edge) const
{
    if (std::find(m_innerEdges.begin(), m_innerEdges.end(), edge) != m_innerEdges.end())
        return true;

    return false;
}


bool shell::Edge::contains(const pmg::Vert* vert) const
{
    if (verts[0]->attachedVert == vert ||
        verts[1]->attachedVert == vert ||
        std::find(m_innerVerts.begin(), m_innerVerts.end(), vert) != m_innerVerts.end())
        return true;

    return false;
}




shell::Edge::Edge(const shell::Vert* vert0, const shell::Vert* vert1)
{
    verts[0] = const_cast<shell::Vert*>(vert0);
    verts[1] = const_cast<shell::Vert*>(vert1);
}
