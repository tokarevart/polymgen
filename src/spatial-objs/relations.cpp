#include "relations.h"


pmg::Edge* pmg::relations::adjacentByEdge(const pmg::Face* face0, const pmg::Face* face1)
{
    size_t inters = 0;
    pmg::Edge* res = nullptr;
    for (auto& edge0 : face0->edges)
        for (auto& edge1 : face1->edges)
            if (edge0 == edge1)
            {
                inters++;
                res = edge0;
                break;
            }

    return inters == 1 ? res : nullptr;
}


bool pmg::relations::contains(const pmg::front::Face* fFace, const pmg::front::Edge* fEdge)
{
    if (fFace->fEdges[0] == fEdge ||
        fFace->fEdges[1] == fEdge ||
        fFace->fEdges[2] == fEdge)
        return true;

    return false;
}


bool pmg::relations::contains(const front::Face* fFace, const pmg::Vert* vert)
{
    return fFace->face->contains(vert);
}


bool pmg::relations::contains(const front::Face* fFace, const pmg::Edge* edge)
{
    return fFace->face->contains(edge);
}
