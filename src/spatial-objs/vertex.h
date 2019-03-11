#pragma once
#include <stddef.h>
#include <memory>
#include "spatial-objs/shell/shell-facet.h"
#include "spatial-objs/shell/shell-edge.h"
#include "spatial-objs/shell/shell-vertex.h"
#include "helpers/spatial-algs/vec.h"

#include "definitions.h"


namespace pmg {

class Vertex
{
    typedef tva::Vec Vec;
    typedef tva::Point Point;

public:
    size_t globalNum;

    shell::Facet*  belongsToShellFacet  = nullptr;
    shell::Edge*   belongsToShellEdge   = nullptr;
    shell::Vertex* belongsToShellVertex = nullptr;

    const Point& pos() const;
          void   setPos( const Point& newPos );
          void   setPos( double coor0, double coor1, double coor2 );

          double& operator[]( short axis );
    const double& operator[]( short axis ) const;
    Vec operator-( const Vertex& other )        const;
    Vec operator-( const shell::Vertex& other ) const;
    Vertex& operator+=( const Vec& other );
    Vertex& operator-=( const Vec& other );

    Vertex();
    Vertex( double coor0, double coor1, double coor2 );
    Vertex( const Point& position );


private:
    std::unique_ptr<Point> m_pos;
};

} // namespace pmg
