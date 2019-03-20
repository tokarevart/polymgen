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
    shell::Vertex* verts[2];

    const std::vector<pmg::Edge*>&   innerEdges() const;
    const std::vector<pmg::Vertex*>& innerVerts() const;

    void segmentize(real_t preferredLen);

    real_t magnitude()    const;
    real_t sqrMagnitude() const;

    bool contains(const shell::Vertex* sVert) const;
    bool contains(const   pmg::Edge*    edge) const;
    bool contains(const   pmg::Vertex*  vert) const;

    Edge(const shell::Vertex* vert0, const shell::Vertex* vert1);


private:
    std::vector<pmg::Edge*>   m_innerEdges;
    std::vector<pmg::Vertex*> m_innerVerts;
};

} // namespace shell
} // namespace pmg
