// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <vector>
#include <memory>
#include "spatial-objs/vertex.h"
#include "helpers/spatial-algs/vec.h"

#include "definitions.h"


namespace pmg {
namespace shell {

class Vertex
{
    using Vec   = tva::Vec;
    using Point = tva::Point;

public:
    pmg::Vertex* attachedVert = nullptr;

    const Point& pos() const;

          double& operator[](unsigned axis);
    const double& operator[](unsigned axis)   const;
    Vec operator-(const shell::Vertex& other) const;
    Vec operator-(const   pmg::Vertex& other) const;

    Vertex();
    Vertex(double coor0, double coor1, double coor2);
    Vertex(const Point& position);


private:
    std::unique_ptr<Point> m_pos;
};

} // namespace shell
} // namespace pmg
