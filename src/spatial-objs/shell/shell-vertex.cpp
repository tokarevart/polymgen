// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/shell/shell-vertex.h"

using namespace pmg;




Vec shell::Vert::operator-(const shell::Vert& other) const
{
    return *m_pos - *other.m_pos;
}


const Vec& shell::Vert::pos() const
{
    return *m_pos;
}


real_t& shell::Vert::operator[](unsigned axis)
{
    return m_pos->coors[axis];
}


const real_t& shell::Vert::operator[](unsigned axis) const
{
    return m_pos->coors[axis];
}


Vec shell::Vert::operator-(const pmg::Vert& other) const
{
    return *m_pos - other.pos();
}




shell::Vert::Vert()
{
    m_pos.reset(new Vec());
}


shell::Vert::Vert(real_t coor0, real_t coor1, real_t coor2)
{
    m_pos.reset(new Vec(coor0, coor1, coor2));
}


shell::Vert::Vert(const Vec& position)
{
    m_pos.reset(new Vec(position));
}
