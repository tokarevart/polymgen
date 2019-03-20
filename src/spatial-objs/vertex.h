// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <stddef.h>
#include <memory>
#include "spatial-objs/shell/shell-face.h"
#include "spatial-objs/shell/shell-edge.h"
#include "spatial-objs/shell/shell-vertex.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class Vertex
{
public:
    size_t globalNum;

    shell::Face*  belongsToShellFace  = nullptr;
    shell::Edge*   belongsToShellEdge   = nullptr;
    shell::Vertex* belongsToShellVertex = nullptr;

    const Vec& pos() const;
          void   setPos( const Vec& newPos );
          void   setPos( real_t coor0, real_t coor1, real_t coor2 );

          real_t& operator[]( short axis );
    const real_t& operator[]( short axis ) const;
    Vec operator-( const Vertex& other )        const;
    Vec operator-( const shell::Vertex& other ) const;
    Vertex& operator+=( const Vec& other );
    Vertex& operator-=( const Vec& other );

    Vertex();
    Vertex( real_t coor0, real_t coor1, real_t coor2 );
    Vertex( const Vec& position );


private:
    std::unique_ptr<Vec> m_pos;
};

} // namespace pmg
