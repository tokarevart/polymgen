// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "spatial-objs/surface/surface-vertex.h"

using namespace pmg;


vec3 surface::Vert::operator-(const surface::Vert& other) const
{
    return *m_pos - *other.m_pos;
}


const vec3& surface::Vert::pos() const
{
    return *m_pos;
}


real_t& surface::Vert::operator[](unsigned axis)
{
    return m_pos->x[axis];
}


const real_t& surface::Vert::operator[](unsigned axis) const
{
    return m_pos->x[axis];
}


vec3 surface::Vert::operator-(const pmg::Vert& other) const
{
    return *m_pos - other.pos();
}




surface::Vert::Vert()
{
    m_pos.reset(new vec3());
}


surface::Vert::Vert(real_t x0, real_t x1, real_t x2)
{
    m_pos.reset(new vec3(x0, x1, x2));
}


surface::Vert::Vert(const vec3& position)
{
    m_pos.reset(new vec3(position));
}
