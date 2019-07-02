// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <memory>
// TODO: remove recursive class dependencies and then remove extra includes
#include "shell/face.h"
#include "shell/edge.h"
#include "shell/vert.h"
#include "../helpers/spatial/vec.h"
#include "../real-type.h"

#include "../definitions.h"


namespace pmg {

class Vert {
    using vec3 = spt::vec<3, real_t>;

public:
    // TODO: use std::unordered_map where i need instead
    std::size_t global_idx;
    real_t min_adj_tetr_vol = std::numeric_limits<real_t>::max();
    real_t max_adj_tetr_vol = std::numeric_limits<real_t>::min();

    // TODO: use std::unordered_map where i need instead
    shell::Face* belongs_to_sface = nullptr;
    shell::Edge* belongs_to_sedge = nullptr;
    shell::Vert* belongs_to_svert = nullptr;

    const vec3& pos() const {
        return *m_pos;
    }
    vec3& pos() {
        return *m_pos;
    }

    Vert() {
        m_pos = std::make_unique<vec3>();
    }
    Vert(const vec3& position) {
        m_pos = std::make_unique<vec3>(position);
    }


private:
    // TODO: move this to public and make normal std::vec3 not ptr after fixing recursive include
    std::unique_ptr<vec3> m_pos;
};

} // namespace pmg
