// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "core/surface/face.h"
#include "helpers/spatial/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg::surface::front {

class Vert
{
    using pair_ee = std::pair<front::Edge*, front::Edge*>;
    using pair_vv = std::pair<pmg::Vert*, pmg::Vert*>;

public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    pmg::Vert* vert;

    void refreshAngleData();

    real_t angle();
    real_t complexity();
    real_t computeAngle();
    real_t computeComplexity();

    pair_ee findAdjEdges() const;
    pair_vv oppVerts()     const;

    Vert(const surface::Face* relatedSurfaceFace, const pmg::Vert* vert);


private:
    // TODO: it's not good to store it
    surface::Face* m_relatedSurfaceFace;

    real_t m_angle;
    real_t m_complexity;

    bool m_needAngleProcessing      = true;
    bool m_needComplexityProcessing = true;
};

} // namespace pmg::surface::front
