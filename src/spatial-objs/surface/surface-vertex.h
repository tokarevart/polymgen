// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <vector>
#include <memory>
#include "spatial-objs/vert.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace surface {

class Vert
{
public:
    pmg::Vert* attachedVert = nullptr;

    const vec3& pos() const;

    real_t&       operator[]( unsigned axis );
    const real_t& operator[]( unsigned axis )  const;
    vec3 operator-( const surface::Vert& other ) const;
    vec3 operator-( const     pmg::Vert& other ) const;

    Vert();
    Vert( real_t x0, real_t x1, real_t x2 );
    Vert( const vec3& position );


private:
    std::unique_ptr<vec3> m_pos;
};

} // namespace surface
} // namespace pmg
