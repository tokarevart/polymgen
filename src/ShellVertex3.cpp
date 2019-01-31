#include "ShellVertex3.h"


tva::Vec3& ShellVertex3::getPos() const
{
    return *m_pos;
}

Vertex3* ShellVertex3::findAttachedVertex(const std::vector<Vertex3*>& freeNodes)
{
    for (auto& vert : freeNodes)
        if (vert->belongsToShellVertex == this)
            return vert;

    return nullptr;
}

double& ShellVertex3::operator[](int axis)
{
    return m_pos->coors[axis];
}

tva::Vec3 ShellVertex3::operator-(const ShellVertex3& right) const
{
    return *m_pos - *right.m_pos;
}

tva::Vec3 ShellVertex3::operator-(const Vertex3& right) const
{
    return *m_pos - right.getPos();
}

ShellVertex3::ShellVertex3()
{
    m_pos.reset(new tva::Vec3());
}

ShellVertex3::ShellVertex3(double coor0, double coor1, double coor2)
{
    m_pos.reset(new tva::Vec3(coor0, coor1, coor2));
}

ShellVertex3::ShellVertex3(const tva::Vec3& position)
{
    m_pos.reset(new tva::Vec3(position));
}

ShellVertex3::~ShellVertex3() {}