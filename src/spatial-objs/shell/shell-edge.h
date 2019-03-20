// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <vector>
#include "spatial-objs/shell/shell-vertex.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace shell {

class Edge
{
public:
    shell::Vert* verts[2];

    const std::vector<pmg::Edge*>&   innerEdges() const;
    const std::vector<pmg::Vert*>& innerVerts() const;

    void segmentize(real_t preferredLen);

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    bool contains(const shell::Vert* sVert) const;
    bool contains(const   pmg::Edge*    edge) const;
    bool contains(const   pmg::Vert*  vert) const;

    Edge(const shell::Vert* vert0, const shell::Vert* vert1);


private:
    std::vector<pmg::Edge*>   m_innerEdges;
    std::vector<pmg::Vert*> m_innerVerts;
};

} // namespace shell
} // namespace pmg
