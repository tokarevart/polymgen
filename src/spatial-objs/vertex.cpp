// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/vertex.h"




const Vec& pmg::Vertex::pos() const
{
    return *m_pos;
}


void pmg::Vertex::setPos(const Vec& newPos)
{
    *m_pos = newPos;
}


void pmg::Vertex::setPos(real_t coor0, real_t coor1, real_t coor2)
{
    m_pos->coors[0] = coor0;
    m_pos->coors[1] = coor1;
    m_pos->coors[2] = coor2;
}




real_t& pmg::Vertex::operator[](short axis)
{
    return m_pos->coors[axis];
}


const real_t& pmg::Vertex::operator[](short axis) const
{
    return m_pos->coors[axis];
}


Vec pmg::Vertex::operator-(const pmg::Vertex& other) const
{
    return *m_pos - *other.m_pos;
}


Vec pmg::Vertex::operator-(const shell::Vertex& other) const
{
    return *m_pos - other.pos();
}


pmg::Vertex& pmg::Vertex::operator+=(const Vec& other)
{
    if (belongsToShellVertex)
    {
        return *this;
    }
    else if (belongsToShellEdge)
    {
        (*m_pos) += Vec(other).project(*belongsToShellEdge->verts[0] - *belongsToShellEdge->verts[1]);
        return *this;
    }
    else if (belongsToShellFacet)
    {
        (*m_pos) += Vec(other).project(
            *belongsToShellFacet->edges[0]->verts[1] - *belongsToShellFacet->edges[0]->verts[0],
            *belongsToShellFacet->edges[1]->verts[1] - *belongsToShellFacet->edges[1]->verts[0]);
        return *this;
    }
    else
    {
        return *this;
    }
}


pmg::Vertex& pmg::Vertex::operator-=(const Vec& other)
{
    if (belongsToShellVertex)
    {
        return *this;
    }
    else if (belongsToShellEdge)
    {
        (*m_pos) -= Vec(other).project(*belongsToShellEdge->verts[0] - *belongsToShellEdge->verts[1]);
        return *this;
    }
    else
    {
        (*m_pos) -= other;
        return *this;
    }
}




pmg::Vertex::Vertex()
{
    m_pos.reset(new Vec());
}


pmg::Vertex::Vertex(real_t coor0, real_t coor1, real_t coor2)
{
    m_pos.reset(new Vec(coor0, coor1, coor2));
}


pmg::Vertex::Vertex(const Vec& position)
{
    m_pos.reset(new Vec(position));
}
