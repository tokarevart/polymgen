#include "relations.h"


pmg::Edge* pmg::relations::adjacentByEdge(const pmg::Face* face0, const pmg::Face* face1) {
    std::size_t inters = 0;
    pmg::Edge* res = nullptr;
    for (auto& edge0 : face0->edges)
        for (auto& edge1 : face1->edges)
            if (edge0 == edge1) {
                inters++;
                res = edge0;
                break;
            }

    return inters == 1 ? res : nullptr;
}
