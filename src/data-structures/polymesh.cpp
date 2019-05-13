// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "data-structures/polymesh.h"

using pmg::PolyMesh;




bool PolyMesh::empty() const
{
    return verts.empty()
        && tetrs.empty()
        && polyhs.empty();
}


void PolyMesh::clear()
{
    verts.clear();
    tetrs.clear();
    polyhs.clear();
}


PolyMesh& PolyMesh::operator=(PolyMesh&& other) noexcept
{
    verts  = std::move(other.verts);
    tetrs  = std::move(other.tetrs);
    polyhs = std::move(other.polyhs);
    return *this;
}


PolyMesh& PolyMesh::operator=(const PolyMesh& other)
{
    verts  = other.verts;
    tetrs  = other.tetrs;
    polyhs = other.polyhs;
    return *this;
}




PolyMesh::PolyMesh(PolyMesh&& other) noexcept
    : verts(std::move(other.verts)), tetrs(std::move(other.tetrs)), polyhs(std::move(other.polyhs)) {}


PolyMesh::PolyMesh(const PolyMesh& other)
{
    verts  = other.verts;
    tetrs  = other.tetrs;
    polyhs = other.polyhs;
}


PolyMesh::PolyMesh() {}
