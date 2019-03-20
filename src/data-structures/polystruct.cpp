// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "data-structures/polystruct.h"

using namespace polygen;




real_t& PolyStruct::VertPos::operator[](unsigned i)
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return z;
    }
}


real_t PolyStruct::VertPos::operator[](unsigned i) const
{
    switch (i)
    {
    case 0: return x;
    case 1: return y;
    case 2: return z;
    default: return z;
    }
}




size_t& PolyStruct::Face::operator[](unsigned i)
{
    switch (i)
    {
    case 0: return vert0;
    case 1: return vert1;
    case 2: return vert2;
    default: return vert2;
    }
}


size_t PolyStruct::Face::operator[](unsigned i) const
{
    switch (i)
    {
    case 0: return vert0;
    case 1: return vert1;
    case 2: return vert2;
    default: return vert2;
    }
}




void PolyStruct::PolyStruct::clear()
{
    verts.clear();
    faces.clear();
    polyhedrons.clear();
}




PolyStruct& PolyStruct::PolyStruct::operator=(PolyStruct&& other) noexcept
{
    verts  = std::move(other.verts);
    faces  = std::move(other.faces);
    polyhedrons = std::move(other.polyhedrons);
    return *this;
}


PolyStruct& PolyStruct::PolyStruct::operator=(const PolyStruct& other)
{
    verts       = other.verts;
    faces       = other.faces;
    polyhedrons = other.polyhedrons;
    return *this;
}




PolyStruct::PolyStruct::PolyStruct(PolyStruct&& other) noexcept
    : verts(std::move(other.verts)), faces(std::move(other.faces)), polyhedrons(std::move(other.polyhedrons)) {}


PolyStruct::PolyStruct::PolyStruct(const PolyStruct& other)
{
    verts       = other.verts;
    faces       = other.faces;
    polyhedrons = other.polyhedrons;
}


PolyStruct::PolyStruct::PolyStruct() {}
