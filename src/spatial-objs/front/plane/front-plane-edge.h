// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#pragma once
#include "spatial-objs/shell/shell-facet.h"
#include "helpers/spatial-algs/vec.h"

#include "definitions.h"


namespace pmg {
namespace front {
namespace plane {

class Edge
{
    typedef tva::Vec   Vec;
    typedef tva::Point Point;

public:
    pmg::Edge* edge;
    Vec normal;

    Vec   computeNormal();
    Point computeCenter();

    static bool isAdj( const Edge* edge0, const Edge* edge1 );

    Edge( const shell::Facet* relatedShellFacet, const pmg::Edge* edge );


private:
    shell::Facet* m_relatedShellFacet;
};

} // namespace plane
} // namespace front
} // namespace pmg
