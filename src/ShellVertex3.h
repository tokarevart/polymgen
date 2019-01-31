#pragma once
#include <vector>
#include <memory>
#include "Vertex3.h"
#include "helpers/spatialalgs/Vec3.h"

class Vertex3;
namespace tva { struct tva::Vec3; }


class ShellVertex3
{
public:
    tva::Vec3& getPos() const;

    Vertex3* findAttachedVertex(
        const std::vector<Vertex3*>& freeNodes);

    double& operator[](int axis);
    tva::Vec3 operator-(const ShellVertex3& right) const;
    tva::Vec3 operator-(const Vertex3& right) const;

    ShellVertex3();
    ShellVertex3(double coor0, double coor1, double coor2);
    ShellVertex3(const tva::Vec3& position);
    ~ShellVertex3();


private:
    std::unique_ptr<tva::Vec3> m_pos;
};