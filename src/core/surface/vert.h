// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <vector>
#include <memory>
#include "core/vert.h"
#include "helpers/spatial/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg::surface {

class Vert
{
public:
    pmg::Vert* attachedVert = nullptr;

    const spt::vec3& pos() const;

    real_t&       operator[]( unsigned axis );
    const real_t& operator[]( unsigned axis )  const;
    spt::vec3 operator-( const surface::Vert& other ) const;
    spt::vec3 operator-( const     pmg::Vert& other ) const;

    Vert();
    Vert( real_t x0, real_t x1, real_t x2 );
    Vert( const spt::vec3& position );


private:
    std::unique_ptr<spt::vec3> m_pos;
};

} // namespace pmg::surface
