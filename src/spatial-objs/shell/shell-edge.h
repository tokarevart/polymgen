#pragma once
#include <vector>
#include "spatial-objs/shell/shell-vertex.h"

#include "definitions.h"


namespace pmg {
namespace shell {

class Edge
{
public:
    shell::Vertex* verts[2];

    const std::vector<pmg::Edge*>&   innerEdges() const;
    const std::vector<pmg::Vertex*>& innerVerts() const;

    void segmentize(double preferredLen);

    double magnitude()    const;
    double sqrMagnitude() const;

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
