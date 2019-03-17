// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/shell/shell-vertex.h"

using namespace pmg;
using tva::Vec;
using tva::Point;




Vec shell::Vertex::operator-(const shell::Vertex& other) const
{
    return *m_pos - *other.m_pos;
}


const Point& shell::Vertex::pos() const
{
    return *m_pos;
}


double& shell::Vertex::operator[](unsigned axis)
{
    return m_pos->coors[axis];
}


const double& shell::Vertex::operator[](unsigned axis) const
{
    return m_pos->coors[axis];
}


Vec shell::Vertex::operator-(const pmg::Vertex& other) const
{
    return *m_pos - other.pos();
}




shell::Vertex::Vertex()
{
    m_pos.reset(new Point());
}


shell::Vertex::Vertex(double coor0, double coor1, double coor2)
{
    m_pos.reset(new Vec(coor0, coor1, coor2));
}


shell::Vertex::Vertex(const Point& position)
{
    m_pos.reset(new Point(position));
}
