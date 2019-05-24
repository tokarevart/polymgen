// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "data-structures/polyshell.h"

using psg::PolyShell;




bool PolyShell::empty() const
{
    return verts.empty()
        && faces.empty()
        && polyhs.empty();
}


void PolyShell::clear()
{
    verts.clear();
    faces.clear();
    polyhs.clear();
}


PolyShell& PolyShell::operator=(PolyShell&& other) noexcept
{
    verts  = std::move(other.verts);
    faces  = std::move(other.faces);
    polyhs = std::move(other.polyhs);
    return *this;
}


PolyShell& PolyShell::operator=(const PolyShell& other)
{
    verts  = other.verts;
    faces  = other.faces;
    polyhs = other.polyhs;
    return *this;
}




PolyShell::PolyShell(PolyShell&& other) noexcept
    : verts(std::move(other.verts)), faces(std::move(other.faces)), polyhs(std::move(other.polyhs)) {}


PolyShell::PolyShell(const PolyShell& other)
{
    verts  = other.verts;
    faces  = other.faces;
    polyhs = other.polyhs;
}


PolyShell::PolyShell() {}
