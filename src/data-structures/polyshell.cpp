// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "data-structures/polyshell.h"

using namespace psg;




real_t& PolyShell::VertPos::operator[](unsigned i)
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return z;
    }
}


real_t PolyShell::VertPos::operator[](unsigned i) const
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return z;
    }
}




size_t& PolyShell::Face::operator[](unsigned i)
{
    switch (i)
    {
    case 0: return vert0;
    case 1: return vert1;
    case 2: return vert2;
    default: return vert2;
    }
}


size_t PolyShell::Face::operator[](unsigned i) const
{
    switch (i)
    {
    case 0: return vert0;
    case 1: return vert1;
    case 2: return vert2;
    default: return vert2;
    }
}




void PolyShell::PolyShell::clear()
{
    verts.clear();
    faces.clear();
    polyhedrons.clear();
}




PolyShell& PolyShell::PolyShell::operator=(PolyShell&& other) noexcept
{
    verts  = std::move(other.verts);
    faces  = std::move(other.faces);
    polyhedrons = std::move(other.polyhedrons);
    return *this;
}


PolyShell& PolyShell::PolyShell::operator=(const PolyShell& other)
{
    verts       = other.verts;
    faces       = other.faces;
    polyhedrons = other.polyhedrons;
    return *this;
}




PolyShell::PolyShell::PolyShell(PolyShell&& other) noexcept
    : verts(std::move(other.verts)), faces(std::move(other.faces)), polyhedrons(std::move(other.polyhedrons)) {}


PolyShell::PolyShell::PolyShell(const PolyShell& other)
{
    verts       = other.verts;
    faces       = other.faces;
    polyhedrons = other.polyhedrons;
}


PolyShell::PolyShell::PolyShell() {}
