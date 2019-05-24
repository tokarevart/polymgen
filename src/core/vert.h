// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstddef>
#include <memory>
// TODO: remove recursive class dependencies and then remove extra includes
#include "core/shell/face.h"
#include "core/shell/edge.h"
#include "core/shell/vert.h"
#include "helpers/spatial/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Vert
{
public:
    // TODO: use std::map where i need instead
    std::size_t globalIdx;
    real_t minAdjTetrVol = std::numeric_limits<real_t>::max();
    real_t maxAdjTetrVol = std::numeric_limits<real_t>::min();

    // TODO: use std::map where i need instead
    shell::Face* belongsToSFace = nullptr;
    shell::Edge* belongsToSEdge = nullptr;
    shell::Vert* belongsToSVert = nullptr;

    const spt::vec3& pos() const;
    void  setPos( const spt::vec3& newPos );
    void  setPos( real_t coor0, real_t coor1, real_t coor2 );

    real_t&       operator[]( unsigned axis );
    const real_t& operator[]( unsigned axis ) const;
    spt::vec3 operator-( const pmg::Vert&   other ) const;
    spt::vec3 operator-( const shell::Vert& other ) const;
    Vert& operator+=( const spt::vec3& other );
    Vert& operator-=( const spt::vec3& other );

    Vert();
    Vert( real_t coor0, real_t coor1, real_t coor2 );
    Vert( const spt::vec3& position );


private:
    std::unique_ptr<spt::vec3> m_pos;
};

} // namespace pmg
