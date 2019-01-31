#pragma once
#include <memory>
#include "ShellFacet3.h"
#include "ShellEdge3.h"
#include "ShellVertex3.h"
#include "Helpers/SpatialAlgs/Vec3.h"

class ShellFacet3;
class ShellEdge3;
class ShellVertex3;
namespace tva { struct Vec3; }


class Vertex3
{
public:
    size_t globalNum;

    ShellFacet3*  belongsToShellFacet  = nullptr;
    ShellEdge3*   belongsToShellEdge   = nullptr;
    ShellVertex3* belongsToShellVertex = nullptr;
    
    const tva::Point3& getPos() const;
          void         setPos(const tva::Point3& newPos);
          void         setPos(double coor0, double coor1, double coor2);
          
    double& operator[](int axis);
    tva::Vec3 operator-(const Vertex3& right) const;
    tva::Vec3 operator-(const ShellVertex3& right) const;
    Vertex3& operator+=(const tva::Vec3& right);
    Vertex3& operator-=(const tva::Vec3& right);

    Vertex3();
    Vertex3(double coor0, double coor1, double coor2);
    Vertex3(const tva::Point3& position);
    ~Vertex3();


private:
    std::unique_ptr<tva::Vec3> m_pos;
};