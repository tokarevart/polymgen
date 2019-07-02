// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include "../face.h"
#include "../../../helpers/spatial/vec.h"
#include "../../../real-type.h"

#include "../../../definitions.h"


namespace pmg::surface::front {

class Vert {
    using pair_ee = std::pair<front::Edge*, front::Edge*>;
    using pair_vv = std::pair<pmg::Vert*, pmg::Vert*>;

public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    pmg::Vert* x;

    void refresh_angle_data();

    real_t angle();
    real_t complexity();
    real_t compute_angle();
    real_t compute_complexity();

    pair_ee adj_edges() const;
    pair_vv opp_verts() const;

    Vert(const surface::Face* relatedSurfaceFace, const pmg::Vert* vert);


private:
    // TODO: it's not good to store it
    surface::Face* m_related_surface_face;

    real_t m_angle = static_cast<real_t>(0.0);
    real_t m_complexity = static_cast<real_t>(0.0);

    bool m_need_angle_processing = true;
    bool m_need_complexity_processing = true;
};

} // namespace pmg::surface::front
