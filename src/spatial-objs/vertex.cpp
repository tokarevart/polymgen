#include "spatial-objs/vertex.h"

using tva::Vec;
using tva::Point;




const Point& pmg::Vertex::pos() const
{
    return *m_pos;
}


void pmg::Vertex::setPos(const Point& newPos)
{
    *m_pos = newPos;
}


void pmg::Vertex::setPos(double coor0, double coor1, double coor2)
{
    m_pos->coors[0] = coor0;
    m_pos->coors[1] = coor1;
    m_pos->coors[2] = coor2;
}




double& pmg::Vertex::operator[](short axis)
{
    return m_pos->coors[axis];
}


const double& pmg::Vertex::operator[](short axis) const
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
    m_pos.reset(new Point());
}


pmg::Vertex::Vertex(double coor0, double coor1, double coor2)
{
    m_pos.reset(new Point(coor0, coor1, coor2));
}


pmg::Vertex::Vertex(const Point& position)
{
    m_pos.reset(new Point(position));
}
