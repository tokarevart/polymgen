// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/shell/shell-facet.h"
#include "helpers/spatial-algs/vec.h"

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

    double complexity();
    double angleExCos();
    double computeComplexity();
    double computeAngleExCos();
    double computeAngle();

    pair_ee findAdjEdges() const;
    pair_vv findOppVerts() const;

    Vertex(const shell::Facet* relatedShellFacet, const pmg::Vertex* vert);


private:
    shell::Facet* m_relatedShellFacet;

    double m_exCos;
    double m_complexity;

    bool m_needExCosProcessing      = true;
    bool m_needComplexityProcessing = true;
};

} // namespace plane
} // namespace front
} // namespace pmg
