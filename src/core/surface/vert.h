// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <vector>
#include <memory>
#include "../vert.h"
#include "../../helpers/spatial/vec.h"
#include "../../real-type.h"

#include "../../definitions.h"


namespace pmg::surface {

class Vert
{
public:
    pmg::Vert* attachedVert = nullptr;

    const spt::vec3& pos() const
    {
        return *m_pos;
    }

    Vert()
    {
        m_pos = std::make_unique<spt::vec3>();
    }
    Vert( const spt::vec3& position )
    {
        m_pos = std::make_unique<spt::vec3>(position);
    }


private:
    // TODO: move this to public and make normal std::vec3 not ptr after fixing recursive include
    std::unique_ptr<spt::vec3> m_pos;
};

} // namespace pmg::surface
