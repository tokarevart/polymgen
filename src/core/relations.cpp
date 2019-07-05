#include "relations.h"

using namespace pmg;


front::Edge* adjacent_by_fedge(const front::Face* fface0, const front::Face* fface1) {
    std::size_t inters = 0;
    front::Edge* res = nullptr;
    for (auto& fedge0 : fface0->front_edges)
        for (auto& fedge1 : fface1->front_edges)
            if (fedge0 == fedge1 && fedge0 && fedge1) {
                if (++inters > 1)
                    return nullptr;
                res = fedge0;
                break;
            }

    return res;
}


pmg::Edge* relations::adjacent_by_edge(const pmg::Face* face0, const pmg::Face* face1) {
    std::size_t inters = 0;
    pmg::Edge* res = nullptr;
    for (auto& edge0 : face0->edges)
        for (auto& edge1 : face1->edges)
            if (edge0 == edge1) {
                if (++inters > 1)
                    return nullptr;
                res = edge0;
                break;
            }

    return res;
}
