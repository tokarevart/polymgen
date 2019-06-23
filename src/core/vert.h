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
public:
    // TODO: use std::unordered_map where i need instead
    std::size_t globalIdx;
    real_t minAdjTetrVol = std::numeric_limits<real_t>::max();
    real_t maxAdjTetrVol = std::numeric_limits<real_t>::min();

    // TODO: use std::unordered_map where i need instead
    shell::Face* belongsToSFace = nullptr;
    shell::Edge* belongsToSEdge = nullptr;
    shell::Vert* belongsToSVert = nullptr;

    const spt::vec3& pos() const {
        return *m_pos;
    }
    spt::vec3& pos() {
        return *m_pos;
    }

    Vert() {
        m_pos = std::make_unique<spt::vec3>();
    }
    Vert(const spt::vec3& position) {
        m_pos = std::make_unique<spt::vec3>(position);
    }


private:
    // TODO: move this to public and make normal std::vec3 not ptr after fixing recursive include
    std::unique_ptr<spt::vec3> m_pos;
};

} // namespace pmg
