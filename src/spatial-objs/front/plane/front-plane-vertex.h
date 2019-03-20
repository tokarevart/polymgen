// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/shell/shell-face.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace front {
namespace plane {

class Vertex
{
    using FrPlEdge   = front::plane::Edge;
    using FrPlVertex = front::plane::Vertex;
    using pair_ee  = std::pair<FrPlEdge*, FrPlEdge*>;
    using pair_vv  = std::pair<pmg::Vertex*, pmg::Vertex*>;

public:
    pmg::Vertex* vert;

    void refreshAngleData();

    real_t complexity();
    real_t angleExCos();
    real_t computeComplexity();
    real_t computeAngleExCos();
    real_t computeAngle();

    pair_ee findAdjEdges() const;
    pair_vv findOppVerts() const;

    Vertex(const shell::Face* relatedShellFace, const pmg::Vertex* vert);


private:
    shell::Face* m_relatedShellFace;

    real_t m_exCos;
    real_t m_complexity;

    bool m_needExCosProcessing      = true;
    bool m_needComplexityProcessing = true;
};

} // namespace plane
} // namespace front
} // namespace pmg
