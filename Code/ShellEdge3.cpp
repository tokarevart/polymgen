#include "ShellEdge3.h"


const double ShellEdge3::magnitude() const
{
    return sqrt(sqrMagnitude());
}

const double ShellEdge3::sqrMagnitude() const
{
    Vec3 buf = *vertexes[1] - *vertexes[0];
    return Vec3::dotProduct(buf, buf);
}

const bool ShellEdge3::contains(const ShellVertex3& vertex) const
{
    if (vertexes[0] == &vertex ||
        vertexes[1] == &vertex)
        return true;

    return false;
}

ShellEdge3::ShellEdge3()
{
    vertexes[0] = nullptr;
    vertexes[1] = nullptr;
}

ShellEdge3::ShellEdge3(
    ShellVertex3& vertex0, 
    ShellVertex3& vertex1)
{
    vertexes[0] = &vertex0;
    vertexes[1] = &vertex1;
}

ShellEdge3::~ShellEdge3() {}