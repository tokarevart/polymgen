#pragma once
#include "ShellVertex3.h"

class ShellVertex3;


class ShellEdge3
{
public:
    ShellVertex3* verts[2];

    double magnitude() const;
    double sqrMagnitude() const;

    bool contains(const ShellVertex3* vert) const;

    ShellEdge3();
    ShellEdge3(const ShellVertex3* vert0, const ShellVertex3* vert1);
    ~ShellEdge3();
};