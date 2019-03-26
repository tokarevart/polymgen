// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/shell/shell-face.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace shell {
namespace front {

class Vert
{
    using pair_ee  = std::pair<front::Edge*, front::Edge*>;
    using pair_vv  = std::pair<pmg::Vert*, pmg::Vert*>;

public:
    pmg::Vert* vert;

    void refreshAngleData();

    real_t complexity();
    real_t angleExCos();
    real_t computeComplexity();
    real_t computeAngleExCos();
    real_t computeAngle();

    pair_ee findAdjEdges() const;
    pair_vv oppVerts() const;

    Vert(const shell::Face* relatedShellFace, const pmg::Vert* vert);


private:
    shell::Face* m_relatedShellFace;

    real_t m_exCos;
    real_t m_complexity;

    bool m_needExCosProcessing      = true;
    bool m_needComplexityProcessing = true;
};

} // namespace front
} // namespace shell
} // namespace pmg
