#include "Vertex3.h"


const tva::Point3& Vertex3::getPos() const
{
    return *m_pos;
}


void Vertex3::setPos(const tva::Point3& newPos)
{
    *m_pos = newPos;
}


void Vertex3::setPos(double coor0, double coor1, double coor2)
{
    m_pos->coors[0] = coor0;
    m_pos->coors[1] = coor1;
    m_pos->coors[2] = coor2;
}




double& Vertex3::operator[](int axis)
{
    return m_pos->coors[axis];
}


tva::Vec3 Vertex3::operator-(const Vertex3& right) const
{
    return *m_pos - *right.m_pos;
}


tva::Vec3 Vertex3::operator-(const ShellVertex3& right) const
{
    return *m_pos - right.getPos();
}


Vertex3& Vertex3::operator+=(const tva::Vec3& right)
{
    if (belongsToShellVertex)
    {
        return *this;
    }
    else if (belongsToShellEdge)
    {
        (*m_pos) += tva::Vec3(right).project(*belongsToShellEdge->verts[0] - *belongsToShellEdge->verts[1]);
        return *this;
    }
    else if (belongsToShellFacet)
    {
        (*m_pos) += tva::Vec3(right).project(
            *belongsToShellFacet->edges[0]->verts[1] - *belongsToShellFacet->edges[0]->verts[0],
            *belongsToShellFacet->edges[1]->verts[1] - *belongsToShellFacet->edges[1]->verts[0]);
        return *this;
    }
    else
    {
        return *this;
    }
}


Vertex3& Vertex3::operator-=(const tva::Vec3& right)
{
    if (belongsToShellVertex)
    {
        return *this;
    }
    else if (belongsToShellEdge)
    {
        (*m_pos) -= tva::Vec3(right).project(*belongsToShellEdge->verts[0] - *belongsToShellEdge->verts[1]);
        return *this;
    }
    else
    {
        (*m_pos) -= right;
        return *this;
    }
}




Vertex3::Vertex3()
{
    m_pos.reset(new tva::Point3());
}


Vertex3::Vertex3(double coor0, double coor1, double coor2)
{
    m_pos.reset(new tva::Point3(coor0, coor1, coor2));
}


Vertex3::Vertex3(const tva::Point3& position)
{
    m_pos.reset(new tva::Point3(position));
}


Vertex3::~Vertex3() {}