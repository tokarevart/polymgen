// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <memory>
#include "spatial-objs/shell/shell-face.h"
#include "spatial-objs/shell/shell-edge.h"
#include "spatial-objs/shell/shell-vertex.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Vert
{
public:
    std::size_t globalIdx;
    real_t minAdjTetrVol = std::numeric_limits<real_t>::max();
    real_t maxAdjTetrVol = std::numeric_limits<real_t>::min();

    shell::Face* belongsToSFace = nullptr;
    shell::Edge* belongsToSEdge = nullptr;
    shell::Vert* belongsToSVert = nullptr;

    const vec3& pos() const;
    void  setPos( const vec3& newPos );
    void  setPos( real_t coor0, real_t coor1, real_t coor2 );

    real_t&       operator[]( unsigned axis );
    const real_t& operator[]( unsigned axis ) const;
    vec3 operator-( const pmg::Vert&   other ) const;
    vec3 operator-( const shell::Vert& other ) const;
    Vert& operator+=( const vec3& other );
    Vert& operator-=( const vec3& other );

    Vert();
    Vert( real_t coor0, real_t coor1, real_t coor2 );
    Vert( const vec3& position );


private:
    std::unique_ptr<vec3> m_pos;
};

} // namespace pmg
