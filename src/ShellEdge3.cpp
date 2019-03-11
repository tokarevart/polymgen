#include "ShellEdge3.h"
#include <cmath>


double ShellEdge3::magnitude() const
{
    return std::sqrt(sqrMagnitude());
}

double ShellEdge3::sqrMagnitude() const
{
    tva::Vec3 buf = *verts[1] - *verts[0];
    return tva::Vec3::dotProduct(buf, buf);
}

bool ShellEdge3::contains(const ShellVertex3* vert) const
{
    if (verts[0] == vert ||
        verts[1] == vert)
        return true;

    return false;
}

ShellEdge3::ShellEdge3()
{
    verts[0] = nullptr;
    verts[1] = nullptr;
}

ShellEdge3::ShellEdge3(const ShellVertex3* vert0, const ShellVertex3* vert1)
{
    verts[0] = const_cast<ShellVertex3*>(vert0);
    verts[1] = const_cast<ShellVertex3*>(vert1);
}

ShellEdge3::~ShellEdge3() {}
