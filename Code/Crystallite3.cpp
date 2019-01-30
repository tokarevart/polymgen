#include "Crystallite3.h"
#include <algorithm>
#include <iostream>

#define DET(a, b, c, d) \
        (a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0[0], corner1[0], point[0]) && \
         BETWEEN(corner0[1], corner1[1], point[1]))

#define PI       3.141592653589793
#define PI_DIV_2 1.5707963267948966

#define DEG_70_IN_RADIANS  1.2217304763960307
#define DEG_80_IN_RADIANS  1.3962634015954636
#define DEG_90_IN_RADIANS  PI_DIV_2
#define DEG_100_IN_RADIANS 1.7453292519943295
#define DEG_150_IN_RADIANS 2.6179938779914943
#define DEG_160_IN_RADIANS 2.7925268031909273
#define DEG_180_IN_RADIANS PI

#define COS_DEG_60   0.5
#define COS_DEG_70   0.342020143325668733
#define COS_DEG_80   0.173648177666930348
#define COS_DEG_90   0.0
#define COS_DEG_100 -0.173648177666930348
#define COS_DEG_120 -0.5
#define COS_DEG_140 -0.766044443118978035
#define COS_DEG_150 -0.866025403784438646
#define COS_DEG_160 -0.939692620785908384
#define COS_DEG_180 -1.0

#define SQRT_2_DIV_3             0.8164965809277260
#define ONE_PLUS_SQRT2_DIV_SQRT3 1.3938468501173517

#define NOT_TOO_CLOSE 1e-1


void Crystallite3::debug()
{
    //unique_ptr<FrontFacet3>* adj_ff[2];
    //for (auto f_edge : frontEdges)
    //{
    //    if (*f_edge)
    //        findAdjacentFrontFacets(f_edge, adj_ff[0], adj_ff[1]);
    //}
    //std::cout << " (" << (size_t)adj_ff[0] * 0 << ") ";

    for (auto f_edge0 : frontEdges)
    {
        if (!*f_edge0)
            continue;

        for (auto f_edge1 : frontEdges)
        {
            if (!*f_edge1 || f_edge1 == f_edge0)
                continue;

            if ((*f_edge0)->edge->vertexes[0] == (*f_edge1)->edge->vertexes[0] &&
                (*f_edge0)->edge->vertexes[1] == (*f_edge1)->edge->vertexes[1] ||
                (*f_edge0)->edge->vertexes[0] == (*f_edge1)->edge->vertexes[1] &&
                (*f_edge0)->edge->vertexes[1] == (*f_edge1)->edge->vertexes[0])
                throw std::exception();
        }
    }
}

template <class T>
void Crystallite3::removePtrsToNullptr(vector<unique_ptr<T>*> &vec)
{
    size_t real_objs_num = std::count_if(vec.begin(), vec.end(), [](unique_ptr<T>*& ptr) { return *ptr ? true : false; });
    vector<unique_ptr<T>*> buf_vec(real_objs_num);

    //size_t firts_thr_nums = real_objs_num / 2ull;
    //size_t second_thr_nums = real_objs_num - firts_thr_nums;
    //#pragma omp parallel num_threads(2)
    //{
    //#pragma omp single
    //{
    size_t index = 0ull;
    for (size_t i = 0ull; index < real_objs_num; i++)
    {
        if (*vec[i])
        {
            buf_vec[index] = vec[i];
            index++;
        }
        else
        {
            delete vec[i];
        }
    }
    //}
    //#pragma omp single
    //{
    //size_t index2 = 0ull;
    //for (size_t i = vec.size() - 1ull; index2 < second_thr_nums; i--)
    //{
    //    if (*vec[i])
    //    {
    //        buf_vec[real_objs_num - 1ull - index2] = vec[i];
    //        index2++;
    //    }
    //}
    //}
    //}

    vec = std::move(buf_vec);
}

const bool Crystallite3::shellContainsVertex(const Vertex3 &vert)
{
    if (vert.belongsToShellFacet)
    {
        for (auto &s_facet : shellFacets)
        {
            if (s_facet == vert.belongsToShellFacet)
                return true;
        }
    }
    else if (vert.belongsToShellEdge)
    {
        for (auto &s_edge : shellEdges)
        {
            if (s_edge == vert.belongsToShellEdge)
                return true;
        }
    }
    else if (vert.belongsToShellVertex)
    {
        for (auto &s_edge : shellEdges)
        {
            if (s_edge->vertexes[0] == vert.belongsToShellVertex ||
                s_edge->vertexes[1] == vert.belongsToShellVertex)
                return true;
        }
    }

    return false;
}

void Crystallite3::setStartFront(
    const vector<unique_ptr<Edge3>*> &edges, 
    const vector<unique_ptr<Facet3>*> &facets)
{
    for (auto &edge : edges)
    {
        if (*edge &&
            shellContainsVertex(**(*edge)->vertexes[0]) &&
            shellContainsVertex(**(*edge)->vertexes[1]))
        {
            frontEdges.push_back((new FrontEdge3(edge->get()))->getPtrToUPtr());
        }
    }

    for (auto &facet : facets)
    {
        if (*facet &&
            shellContainsVertex(**(*(*facet)->edges[0])->vertexes[0]) &&
            shellContainsVertex(**(*(*facet)->edges[0])->vertexes[1]) &&
            shellContainsVertex(**(*facet)->findVertexNotIncludedInEdge(**(*facet)->edges[0])))
        {
            frontFacets.push_back((new FrontFacet3(facet->get()))->getPtrToUPtr());
        }
    }
}

void Crystallite3::computeFrontNormals()
{
    for (auto facet : frontFacets)
    {
        if (!*facet)
            continue;

        (*facet)->computeNormal(frontFacets);
    }
}

ShellEdge3* Crystallite3::findShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const
{
    for (auto &s_edge : shellEdges)
        if ((s_edge->vertexes[0] == v0  &&
             s_edge->vertexes[1] == v1) ||
            (s_edge->vertexes[1] == v0  &&
             s_edge->vertexes[0] == v1))
        {
            return s_edge;
        }

    return nullptr;
}

unique_ptr<FrontFacet3>* Crystallite3::findFrontFacet(const unique_ptr<Facet3>* facet)
{
    for (auto &f_facet : frontFacets)
    {
        if (*f_facet &&
            (*f_facet)->facet == facet->get())
            return f_facet;
    }

    return nullptr;
}

unique_ptr<FrontEdge3>* Crystallite3::findFrontEdge(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1)
{
    for (auto f_edge : frontEdges)
    {
        if (*f_edge &&
            (((*f_edge)->edge->vertexes[0] == v0 &&
              (*f_edge)->edge->vertexes[1] == v1) ||
             ((*f_edge)->edge->vertexes[1] == v0 &&
              (*f_edge)->edge->vertexes[0] == v1)))
            return f_edge;
    }

    return nullptr;
}

unique_ptr<FrontEdge3>* Crystallite3::findFrontEdge(const unique_ptr<Edge3>* edge)
{
    for (auto &f_edge : frontEdges)
    {
        if (*f_edge &&
            (*f_edge)->edge == edge->get())
            return f_edge;
    }

    return nullptr;
}

void Crystallite3::findAdjacentFrontFacets(const unique_ptr<FrontEdge3>* frontEdge, unique_ptr<FrontFacet3>*& out_frontFacet0, unique_ptr<FrontFacet3>*& out_frontFacet1) const
{
    bool not_found_yet = true;
    for (auto f_facet : frontFacets)
    {
        if (*f_facet &&
            (*f_facet)->facet->contains(*(*frontEdge)->edge))
        {
            if (not_found_yet)
            {
                out_frontFacet0 = f_facet;
                not_found_yet = false;
            }
            else
            {
                out_frontFacet1 = f_facet;
                return;
            }
        }
    }

    throw std::logic_error("Something went wrong in Crystallite3::findAdjacentFrontFacets");
}

void Crystallite3::findOppositeVertexes(const unique_ptr<FrontEdge3>* frontEdge, unique_ptr<Vertex3>*& out_vert0, unique_ptr<Vertex3>*& out_vert1) const
{
    unique_ptr<FrontFacet3>* around_f_facets[2];
    findAdjacentFrontFacets(frontEdge, around_f_facets[0], around_f_facets[1]);

    out_vert0 = (*around_f_facets[0])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);
    out_vert1 = (*around_f_facets[1])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);
}

unique_ptr<FrontEdge3>* Crystallite3::findOppositeFrontEdge(const unique_ptr<FrontEdge3>* frontEdge) const 
{
    unique_ptr<Vertex3>* opp_verts[2];
    findOppositeVertexes(frontEdge, opp_verts[0], opp_verts[1]);

    for (auto f_edge : frontEdges)
    {
        if (*f_edge &&
            (((*f_edge)->edge->vertexes[0] == opp_verts[0]  &&
              (*f_edge)->edge->vertexes[1] == opp_verts[1]) ||
             ((*f_edge)->edge->vertexes[1] == opp_verts[0]  &&
              (*f_edge)->edge->vertexes[0] == opp_verts[1])))
            return f_edge;
    }

    return nullptr;
}

unique_ptr<FrontFacet3>* Crystallite3::addFrontFacet3(unique_ptr<Edge3>* edge0, unique_ptr<Edge3>* edge1, unique_ptr<Edge3>* edge2)
{
    FrontFacet3* new_f_facet = new FrontFacet3(new Facet3(**edge0, **edge1, **edge2));

    frontFacets.push_back(new_f_facet->getPtrToUPtr());
    innerFacets.push_back(new_f_facet->facet->getPtrToUPtr());

    return new_f_facet->getPtrToUPtr();
}

unique_ptr<FrontEdge3>* Crystallite3::addFrontEdge3(unique_ptr<Vertex3>* vert0, unique_ptr<Vertex3>* vert1)
{
    FrontEdge3* new_f_edge = new FrontEdge3(new Edge3(**vert0, **vert1));

    frontEdges.push_back(new_f_edge->getPtrToUPtr());
    innerEdges.push_back(new_f_edge->edge->getPtrToUPtr());

    return new_f_edge->getPtrToUPtr();
}

void Crystallite3::addFrontEdge3(unique_ptr<FrontEdge3>* frontEdge)
{
    frontEdges.push_back((*frontEdge)->getPtrToUPtr());
    innerEdges.push_back((*frontEdge)->edge->getPtrToUPtr());
}

const bool Crystallite3::vertexInsideFrontCheck(const Vec3& v)
{
    int inters_num = 0;
    for (auto &f_facet : frontFacets)
    {
        if (!*f_facet)
            continue;
        
        Vec3 dir = 0.3333333333333333 * ((*f_facet)->computeCenter() - 3.0 * v);

        for (auto &f_facetj : frontFacets)
        {
            if (*f_facetj && 
                spatialalgs::isRayIntersectTriangle(
                    v, dir,
                    (*(*(*f_facetj)->facet->edges[0])->vertexes[0])->getPosition(),
                    (*(*(*f_facetj)->facet->edges[0])->vertexes[1])->getPosition(),
                    (*(*f_facetj)->facet->findVertexNotIncludedInEdge(**(*f_facetj)->facet->edges[0]))->getPosition()))
            {
                inters_num++;
            }
        }

        break;
    }

    return inters_num % 2 == 1;
}

const bool Crystallite3::lineSegmentGlobalIntersectionCheck(const Vec3& v0, const Vec3& v1)
{
    for (auto &facet : innerFacets)
    {
        if (*facet &&
            spatialalgs::isSegmentIntersectTriangle(
                v0, v1,
                (*(*(*facet)->edges[0])->vertexes[0])->getPosition(),
                (*(*(*facet)->edges[0])->vertexes[1])->getPosition(),
                (*(*facet)->findVertexNotIncludedInEdge(**(*facet)->edges[0]))->getPosition()))
            return true;
    }

    return false;
}

const bool Crystallite3::segmentFrontIntersectionCheck(const Vec3& v0, const Vec3& v1)
{
    for (auto &f_facet : frontFacets)
    {
        if (*f_facet &&
            spatialalgs::isSegmentIntersectTriangle(
                v0, v1,
                (*(*(*f_facet)->facet->edges[0])->vertexes[0])->getPosition(),
                (*(*(*f_facet)->facet->edges[0])->vertexes[1])->getPosition(),
                (*(*f_facet)->facet->findVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->getPosition()))
            return true;
    }

    return false;
}

const bool Crystallite3::edgeGlobalIntersectionCheck(const unique_ptr<Edge3>* edge)
{
    Vec3 delta = 1e-3 * (**(*edge)->vertexes[1] - **(*edge)->vertexes[0]);
    Vec3 segment[2];
    segment[0] = (*(*edge)->vertexes[0])->getPosition() + delta;
    segment[1] = (*(*edge)->vertexes[1])->getPosition() - delta;

    for (auto &facet : innerFacets)
    {
        if (*facet &&
            !(*facet)->contains(**edge) &&
            spatialalgs::isSegmentIntersectTriangle(
                segment[0], segment[1],
                (*(*(*facet)->edges[0])->vertexes[0])->getPosition(), 
                (*(*(*facet)->edges[0])->vertexes[1])->getPosition(), 
                (*(*facet)->findVertexNotIncludedInEdge(**(*facet)->edges[0]))->getPosition()))
            return true;
    }

    return false;
}

const bool XOR(const bool b0, const bool b1)
{
    return (b0 || b1) && !(b0 && b1);
}

const bool Crystallite3::edgeIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge)
{
    unique_ptr<Vertex3>* opp_verts[2];
    findOppositeVertexes(
        frontEdge,
        opp_verts[0],
        opp_verts[1]);
    Vec3 opp_verts_poses[2];
    opp_verts_poses[0] = (*opp_verts[0])->getPosition();
    opp_verts_poses[1] = (*opp_verts[1])->getPosition();
    Vec3 opp_vert0_to_1 = opp_verts_poses[1] - opp_verts_poses[0];
    Vec3 delta_vec_0th_to_1st = 1e-3 * opp_vert0_to_1;

    if (segmentFrontIntersectionCheck(
        opp_verts_poses[0] + delta_vec_0th_to_1st,
        opp_verts_poses[1] - delta_vec_0th_to_1st))
        return true;

    for (auto &f_edge : frontEdges)
    {
        if (*f_edge &&
            (*f_edge)->edge->contains(**opp_verts[0]) &&
            (*f_edge)->edge->contains(**opp_verts[1]))
            return false;
    }

    for (auto &f_edge : frontEdges)
    {
        if (!*f_edge)
            continue;

        bool contains[2];
        contains[0] = (*f_edge)->edge->contains(**opp_verts[0]);
        contains[1] = (*f_edge)->edge->contains(**opp_verts[1]);

        unique_ptr<Vertex3>* vert_buf;
        if (contains[0])
        {
            if ((*f_edge)->edge->vertexes[0] == opp_verts[0])
                vert_buf = (*f_edge)->edge->vertexes[1];
            else
                vert_buf = (*f_edge)->edge->vertexes[0];

            if ((*vert_buf)->getPosition().distanceToSegment(opp_verts_poses[0], opp_verts_poses[1]) < 4e-3 * _preferredLength)
                return true;
        }
        else if (contains[1])
        {
            if ((*f_edge)->edge->vertexes[0] == opp_verts[1])
                vert_buf = (*f_edge)->edge->vertexes[1];
            else
                vert_buf = (*f_edge)->edge->vertexes[0];

            if ((*vert_buf)->getPosition().distanceToSegment(opp_verts_poses[0], opp_verts_poses[1]) < 4e-3 * _preferredLength)
                return true;
        }
        else
        {
            if (spatialalgs::segmentsDistance(
                    opp_verts_poses[0], opp_verts_poses[1],
                    (*(*f_edge)->edge->vertexes[0])->getPosition(), (*(*f_edge)->edge->vertexes[1])->getPosition()) < 4e-3 * _preferredLength)
                return true;
        }
    }

    return false;
}

const bool Crystallite3::facetIntersectionCheck(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1, const Vec3& v2)
{
    for (auto &f_edge : frontEdges)
    {
        if (!*f_edge)
            continue;

        bool contains[2];
        contains[0] = (*f_edge)->edge->contains(**v0);
        contains[1] = (*f_edge)->edge->contains(**v1);

        if (contains[0] || contains[1])
            continue;

        Vec3 first_to_second_delta = NOT_TOO_CLOSE * ((*(*f_edge)->edge->vertexes[1])->getPosition() - (*(*f_edge)->edge->vertexes[0])->getPosition());
        Vec3 f_edge_verts_poses[2];
        f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition() - first_to_second_delta;
        f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition() + first_to_second_delta;

        if (spatialalgs::isSegmentIntersectTriangle(
            f_edge_verts_poses[0], f_edge_verts_poses[1],
            (*v0)->getPosition(), (*v1)->getPosition(), v2))
            return true;
    }

    return false;
}

//const bool Crystallite3::facetIntersectionCheck(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1, const Vec3& v2)
//{
//    Vec3 verts_poses[3];
//    verts_poses[0] = (*v0)->getPosition();
//    verts_poses[1] = (*v1)->getPosition();
//
//    for (auto &f_edge : frontEdges)
//    {
//        if (!*f_edge)
//            continue;
//
//        bool contains[3];
//        contains[0] = (*f_edge)->edge->contains(**v0);
//        contains[1] = (*f_edge)->edge->contains(**v1);
//
//        if (contains[0] && contains[1])
//            continue;
//
//        Vec3 first_to_second_delta = 1e-3 * ((*(*f_edge)->edge->vertexes[1])->getPosition() - (*(*f_edge)->edge->vertexes[0])->getPosition());
//        Vec3 f_edge_verts_poses[2];
//        if (contains[0])
//        {
//            if ((*f_edge)->edge->vertexes[0] == v0)
//            {
//                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition() + first_to_second_delta;
//                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition();
//            }
//            else
//            {
//                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition();
//                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition() - first_to_second_delta;
//            }
//        }
//        else if (contains[1])
//        {
//            if ((*f_edge)->edge->vertexes[0] == v1)
//            {
//                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition() + first_to_second_delta;
//                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition();
//            }
//            else
//            {
//                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition();
//                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition() - first_to_second_delta;
//            }
//        }
//        else
//        {
//            f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition();
//            f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition();
//        }
//
//        if (spatialalgs::isSegmentIntersectTriangle(
//                f_edge_verts_poses[0], f_edge_verts_poses[1],
//                verts_poses[0], verts_poses[1], v2))
//            return true;
//    }
//
//    return false;
//}

const bool Crystallite3::facetIntersectionCheck(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1, const unique_ptr<Vertex3>* v2)
{
    Vec3 verts_poses[3];
    verts_poses[0] = (*v0)->getPosition();
    verts_poses[1] = (*v1)->getPosition();
    verts_poses[2] = (*v2)->getPosition();

    for (auto &f_edge : frontEdges)
    {
        if (!*f_edge)
            continue;

        bool contains[3];
        contains[0] = (*f_edge)->edge->contains(**v0);
        contains[1] = (*f_edge)->edge->contains(**v1);
        contains[2] = (*f_edge)->edge->contains(**v2);

        if ((contains[0] && contains[1]) ||
            (contains[0] && contains[2]) ||
            (contains[1] && contains[2]))
            continue;

        Vec3 first_to_second_delta = 1e-3 * ((*(*f_edge)->edge->vertexes[1])->getPosition() - (*(*f_edge)->edge->vertexes[0])->getPosition());
        Vec3 f_edge_verts_poses[2];
        if (contains[0])
        {
            if ((*f_edge)->edge->vertexes[0] == v0)
            {
                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition() + first_to_second_delta;
                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition();
            }
            else
            {
                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition();
                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition() - first_to_second_delta;
            }
        }
        else if (contains[1])
        {
            if ((*f_edge)->edge->vertexes[0] == v1)
            {
                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition() + first_to_second_delta;
                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition();
            }
            else
            {
                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition();
                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition() - first_to_second_delta;
            }
        }
        else if (contains[2])
        {
            if ((*f_edge)->edge->vertexes[0] == v2)
            {
                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition() + first_to_second_delta;
                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition();
            }
            else
            {
                f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition();
                f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition() - first_to_second_delta;
            }
        }
        else
        {
            f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->getPosition();
            f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->getPosition();
        }

        if (spatialalgs::isSegmentIntersectTriangle(
                f_edge_verts_poses[0], f_edge_verts_poses[1],
                verts_poses[0],    verts_poses[1], verts_poses[2]))
            return true;
    }

    return false;
}

const bool Crystallite3::facetsIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge)
{
    unique_ptr<Vertex3>* opp_verts[2];
    findOppositeVertexes(frontEdge, opp_verts[0], opp_verts[1]);

    Vec3 opp_verts_poses[2];
    opp_verts_poses[0] = (*opp_verts[0])->getPosition();
    opp_verts_poses[1] = (*opp_verts[1])->getPosition();

    //for (auto f_edge : frontEdges)
    //{
    //    if (*f_edge &&
    //        (((*f_edge)->edge->vertexes[0] == opp_verts[0] &&
    //          (*f_edge)->edge->vertexes[1] == opp_verts[1]) ||
    //         ((*f_edge)->edge->vertexes[1] == opp_verts[0] &&
    //          (*f_edge)->edge->vertexes[0] == opp_verts[1])))
    //        return false;
    //}

    if (facetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[0], opp_verts[1]) ||
        facetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], opp_verts[1]))
        return true;

    return false;
}

const bool Crystallite3::insideSimplex3Check(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& vert)
{
    Vec3 vert_to_p0 = p0 - vert;
    Vec3 vert_to_p1 = p1 - vert;
    Vec3 vert_to_p2 = p2 - vert;
    Vec3 vert_to_p3 = p3 - vert;

    double abs_mixed_prods[5];
    abs_mixed_prods[0] = abs(Vec3::mixedProduct(vert_to_p0, vert_to_p2, vert_to_p3));
    abs_mixed_prods[1] = abs(Vec3::mixedProduct(vert_to_p0, vert_to_p1, vert_to_p2));
    abs_mixed_prods[2] = abs(Vec3::mixedProduct(vert_to_p0, vert_to_p1, vert_to_p3));
    abs_mixed_prods[3] = abs(Vec3::mixedProduct(vert_to_p1, vert_to_p2, vert_to_p3));
    abs_mixed_prods[4] = abs(Vec3::mixedProduct(p1 - p0, p2 - p0, p3 - p0));

    /*if (abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4] * 1.0000001)
        return true;
    else
        return false;*/
    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4] * 1.000001;
}

const bool Crystallite3::anyVertexInsidePotentialSimplex3Check(const unique_ptr<FrontEdge3>* frontEdge)
{
    unique_ptr<Vertex3>* opp_verts[2];
    findOppositeVertexes(frontEdge, opp_verts[0], opp_verts[1]);

    Vec3 points[4];
    points[0] = (*opp_verts[0])->getPosition();
    points[1] = (*opp_verts[1])->getPosition();
    points[2] = (*(*frontEdge)->edge->vertexes[0])->getPosition();
    points[3] = (*(*frontEdge)->edge->vertexes[1])->getPosition();
    for (auto &vert : innerVerts)
    {
        if (*vert &&
            vert != opp_verts[0] &&
            vert != opp_verts[1] &&
            vert != (*frontEdge)->edge->vertexes[0] &&
            vert != (*frontEdge)->edge->vertexes[1] &&
            insideSimplex3Check(points[0], points[1], points[2], points[3], (*vert)->getPosition()))
        {
            return true;
        }
    }

    return false;
}

const bool Crystallite3::frontSplitCheck(const unique_ptr<FrontEdge3>* frontEdge)
{
    unique_ptr<FrontEdge3>* opp_f_edge = findOppositeFrontEdge(frontEdge);
    if (!opp_f_edge)
        return false;

    unique_ptr<Vertex3>* opp_f_edge_opp_verts[2];
    findOppositeVertexes(
        opp_f_edge,
        opp_f_edge_opp_verts[0], 
        opp_f_edge_opp_verts[1]);

    if (opp_f_edge_opp_verts[0] == (*frontEdge)->edge->vertexes[0] ||
        opp_f_edge_opp_verts[0] == (*frontEdge)->edge->vertexes[1] ||
        opp_f_edge_opp_verts[1] == (*frontEdge)->edge->vertexes[0] ||
        opp_f_edge_opp_verts[1] == (*frontEdge)->edge->vertexes[1])
        return false;

    return true;
}

const bool Crystallite3::parallelFacetsCheck(const unique_ptr<FrontEdge3>* frontEdge) const
{
    unique_ptr<FrontFacet3>* adj_f_facets[2];
    findAdjacentFrontFacets(frontEdge, adj_f_facets[0], adj_f_facets[1]);

    unique_ptr<Vertex3>* opp_verts[2];
    opp_verts[0] = (*adj_f_facets[0])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);
    opp_verts[1] = (*adj_f_facets[1])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);

    Vec3 plane0[2];
    plane0[0] = **opp_verts[0] - **(*frontEdge)->edge->vertexes[0];
    plane0[1] = **opp_verts[1] - **(*frontEdge)->edge->vertexes[0];
    Vec3 normal0 = Vec3::crossProduct(plane0[0], plane0[1]).normalize();
    Vec3 plane1[2];
    plane1[0] = **opp_verts[0] - **(*frontEdge)->edge->vertexes[1];
    plane1[1] = **opp_verts[1] - **(*frontEdge)->edge->vertexes[1];
    Vec3 normal1 = Vec3::crossProduct(plane1[0], plane1[1]).normalize();

    for (auto &f_facet : frontFacets)
    {
        unique_ptr<Edge3>* inters_reses[2];
        if (*f_facet &&
            f_facet != adj_f_facets[0] &&
            f_facet != adj_f_facets[1])
        {
            inters_reses[0] = Facet3::intersectAlongAnEdge(*(*f_facet)->facet, *(*adj_f_facets[0])->facet);
            inters_reses[1] = Facet3::intersectAlongAnEdge(*(*f_facet)->facet, *(*adj_f_facets[1])->facet);
            if (!XOR(inters_reses[0], inters_reses[1]))
                continue;

            unique_ptr<Vertex3>* f_facet_to_verts[2];
            f_facet_to_verts[0] = (*(*f_facet)->facet->edges[0])->vertexes[0];
            f_facet_to_verts[1] = (*(*f_facet)->facet->edges[0])->vertexes[1];
            unique_ptr<Vertex3>* f_facet_from_vert = (*f_facet)->facet->findVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]);

            Vec3 f_plane[2];
            f_plane[0] = **f_facet_to_verts[0] - **f_facet_from_vert;
            f_plane[1] = **f_facet_to_verts[1] - **f_facet_from_vert;
            Vec3 f_normal = Vec3::crossProduct(f_plane[0], f_plane[1]).normalize();

            if (abs(abs(Vec3::dotProduct(f_normal, normal0)) - 1.0) < 1e-6 ||
                abs(abs(Vec3::dotProduct(f_normal, normal1)) - 1.0) < 1e-6)
            {
                int i = inters_reses[0] ? 0 : 1;

                Vec3 border_verts[2];
                border_verts[0] = (*(*inters_reses[i])->vertexes[0])->getPosition();
                border_verts[1] = (*(*inters_reses[i])->vertexes[1])->getPosition();

                Vec3 main_facet_3rd_vert;
                if ((*inters_reses[i])->contains(**opp_verts[0]))
                    main_facet_3rd_vert = (*opp_verts[1])->getPosition();
                else
                    main_facet_3rd_vert = (*opp_verts[0])->getPosition();

                Vec3 curr_facet_3rd_vert = (*(*f_facet)->facet->findVertexNotIncludedInEdge(**inters_reses[i]))->getPosition();

                Vec3 main_facet_cross = Vec3::crossProduct(main_facet_3rd_vert - border_verts[0], main_facet_3rd_vert - border_verts[1]);
                Vec3 curr_facet_cross = Vec3::crossProduct(curr_facet_3rd_vert - border_verts[0], curr_facet_3rd_vert - border_verts[1]);
                if (Vec3::dotProduct(main_facet_cross, curr_facet_cross) > 0.0)
                    return true;
            }
        }
    }

    return false;
}

const bool Crystallite3::frontContainsOfOnly1Simplex3OrEmpty() const
{
    size_t f_facets_num = 0ull;
    for (auto &f_facet : frontFacets)
    {
        if (*f_facet)
        {
            f_facets_num++;
            if (f_facets_num > 4ull)
                return false;
        }
    }

    return true;
}

const double Crystallite3::computeMinEdgesLength(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    double sqr_magns[6];
    sqr_magns[0] = (p1 - p0).sqrMagnitude();
    sqr_magns[1] = (p2 - p0).sqrMagnitude();
    sqr_magns[2] = (p3 - p0).sqrMagnitude();
    sqr_magns[3] = (p2 - p1).sqrMagnitude();
    sqr_magns[4] = (p3 - p1).sqrMagnitude();
    sqr_magns[5] = (p3 - p2).sqrMagnitude();
    return sqrt(std::min({ sqr_magns[0], sqr_magns[1], sqr_magns[2], sqr_magns[3], sqr_magns[4], sqr_magns[5] }));
}

const double Crystallite3::computeMaxEdgesLength(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    double sqr_magns[6];
    sqr_magns[0] = (p1 - p0).sqrMagnitude();
    sqr_magns[1] = (p2 - p0).sqrMagnitude();
    sqr_magns[2] = (p3 - p0).sqrMagnitude();
    sqr_magns[3] = (p2 - p1).sqrMagnitude();
    sqr_magns[4] = (p3 - p1).sqrMagnitude();
    sqr_magns[5] = (p3 - p2).sqrMagnitude();
    return sqrt(std::max({ sqr_magns[0], sqr_magns[1], sqr_magns[2], sqr_magns[3], sqr_magns[4], sqr_magns[5] }));
}

std::pair<double, double> Crystallite3::computeMinMaxEdgesLengths(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    double sqr_magns[6];
    sqr_magns[0] = (p1 - p0).sqrMagnitude();
    sqr_magns[1] = (p2 - p0).sqrMagnitude();
    sqr_magns[2] = (p3 - p0).sqrMagnitude();
    sqr_magns[3] = (p2 - p1).sqrMagnitude();
    sqr_magns[4] = (p3 - p1).sqrMagnitude();
    sqr_magns[5] = (p3 - p2).sqrMagnitude();
    auto min_max = std::minmax({ sqr_magns[0], sqr_magns[1], sqr_magns[2], sqr_magns[3], sqr_magns[4], sqr_magns[5] });
    min_max.first = sqrt(min_max.first);
    min_max.second = sqrt(min_max.second);
    return min_max;
}

std::pair<double, double> Crystallite3::computeMinMaxEdgesSqrLengths(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    double sqr_magns[6];
    sqr_magns[0] = (p1 - p0).sqrMagnitude();
    sqr_magns[1] = (p2 - p0).sqrMagnitude();
    sqr_magns[2] = (p3 - p0).sqrMagnitude();
    sqr_magns[3] = (p2 - p1).sqrMagnitude();
    sqr_magns[4] = (p3 - p1).sqrMagnitude();
    sqr_magns[5] = (p3 - p2).sqrMagnitude();
    return std::minmax({ sqr_magns[0], sqr_magns[1], sqr_magns[2], sqr_magns[3], sqr_magns[4], sqr_magns[5] });
}

const double Crystallite3::computeSimplex3SimpleQuality(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return sqrt(sqr_min_max.first / sqr_min_max.second);
}

const double Crystallite3::computeSimplex3SimpleSqrQuality(const Vec3 & p0, const Vec3 & p1, const Vec3 & p2, const Vec3 & p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return sqr_min_max.first / sqr_min_max.second;
}

unique_ptr<FrontEdge3>* Crystallite3::currentFrontEdge(double maxExCos) const
{
    double curr_max_excos = -2.0;
    unique_ptr<FrontEdge3>* curr_max_f_edge = nullptr;
    for (size_t i = 0ull; i < frontEdges.size(); i++)
    {
        if (!*frontEdges[i])
            continue;

        double curr_excos = (*frontEdges[i])->getAngleExCos(this);
        if (curr_excos > curr_max_excos &&
            curr_excos < maxExCos)
        {
            curr_max_excos = curr_excos;
            curr_max_f_edge = frontEdges[i];
        }
    }

    return curr_max_f_edge;
}

const bool Crystallite3::exhaustWithoutNewVertexPriorityPredicate(unique_ptr<FrontEdge3>* currentFrontEdge)
{
    if ((*currentFrontEdge)->getAngleExCos(this) > COS_DEG_70)
        return true;

    unique_ptr<Vertex3>* opp_verts[2];
    findOppositeVertexes(currentFrontEdge, opp_verts[0], opp_verts[1]);
    if (findFrontEdge(opp_verts[0], opp_verts[1]) ||
        ((*currentFrontEdge)->getAngleExCos(this) < COS_DEG_70 &&
         (*currentFrontEdge)->getAngleExCos(this) > COS_DEG_100 &&
         (**opp_verts[1] - **opp_verts[0]).sqrMagnitude() <= _preferredLength * _preferredLength))
        return true;

    return false;
}

const bool Crystallite3::exhaustWithNewVertexPriorityPredicate(unique_ptr<FrontEdge3>* currentFrontEdge)
{
    if ((*currentFrontEdge)->getAngleExCos(this) < COS_DEG_120)
        return true;

    return false;
}

const Crystallite3::ExhaustType Crystallite3::exhaustTypeQualityPriorityCalculation(
    unique_ptr<FrontEdge3>* currentFrontEdge, 
    unique_ptr<FrontFacet3>** out_withNWFrontFacet, Vec3** out_withNWNewVertPos)
{
    if (frontSplitCheck(currentFrontEdge))
        return DONT_EXHAUST;
    
    if (parallelFacetsCheck(currentFrontEdge) ||
        edgeIntersectionCheck(currentFrontEdge) ||
        facetsIntersectionCheck(currentFrontEdge) ||
        anyVertexInsidePotentialSimplex3Check(currentFrontEdge))
    {
        return WITH_NEW_VERTEX;
    }

    unique_ptr<Vertex3>* opp_verts[2];
    findOppositeVertexes(currentFrontEdge, opp_verts[0], opp_verts[1]);
    double without_nv_quality = computeSimplex3SimpleSqrQuality(
        (*(*currentFrontEdge)->edge->vertexes[0])->getPosition(),
        (*(*currentFrontEdge)->edge->vertexes[1])->getPosition(),
        (*opp_verts[0])->getPosition(),
        (*opp_verts[1])->getPosition());

    auto f_facet = chooseFrontFacetForExhaustionWithNewVertex(currentFrontEdge);
    Vec3 new_vert_pos;
    if (!tryComputeNewVertexPosition(f_facet, new_vert_pos))
        return DONT_EXHAUST;

    double with_nv_quality = computeSimplex3SimpleSqrQuality(
        (*(*(*f_facet)->facet->edges[0])->vertexes[0])->getPosition(),
        (*(*(*f_facet)->facet->edges[0])->vertexes[1])->getPosition(),
        (*(*f_facet)->facet->findVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->getPosition(),
        new_vert_pos);

    if (without_nv_quality > with_nv_quality)
        return WITHOUT_NEW_VERTEX;

    *out_withNWFrontFacet = f_facet;
    *out_withNWNewVertPos = new Vec3(new_vert_pos);
    return WITH_NEW_VERTEX;
}

Vec3 Crystallite3::computeNormalInSimplex3(unique_ptr<FrontFacet3>* frontFacet, const Vec3& oppositeVertPos)
{
    Vec3 third_pos = (*(*frontFacet)->facet->findVertexNotIncludedInEdge(**(*frontFacet)->facet->edges[0]))->getPosition();
    Vec3 normal = Vec3::crossProduct(
        (*(*(*frontFacet)->facet->edges[0])->vertexes[0])->getPosition() - third_pos,
        (*(*(*frontFacet)->facet->edges[0])->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, oppositeVertPos - third_pos) > 0.0)
        normal *= -1.0;

    return normal;
}

Vec3 Crystallite3::computeNormalInSimplex3(unique_ptr<FrontFacet3>* frontFacet, unique_ptr<Edge3>* oneOfTheRemainingEdges)
{
    Vec3 opposite = (*frontFacet)->facet->contains(**(*oneOfTheRemainingEdges)->vertexes[0]) ?
        (*(*oneOfTheRemainingEdges)->vertexes[1])->getPosition() :
        (*(*oneOfTheRemainingEdges)->vertexes[0])->getPosition();
    Vec3 third_pos = (*(*frontFacet)->facet->findVertexNotIncludedInEdge(**(*frontFacet)->facet->edges[0]))->getPosition();
    Vec3 normal = Vec3::crossProduct(
        (*(*(*frontFacet)->facet->edges[0])->vertexes[0])->getPosition() - third_pos,
        (*(*(*frontFacet)->facet->edges[0])->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - third_pos) > 0.0)
        normal *= -1.0;

    return normal;
}

void Crystallite3::exhaustWithoutNewVertexOppositeEdgeExists(unique_ptr<FrontEdge3>* frontEdge, unique_ptr<FrontEdge3>* oppositeEdge)
{
    unique_ptr<FrontEdge3>* opp_f_edge = oppositeEdge;

    unique_ptr<FrontFacet3>* opp_f_facets[2];
    findAdjacentFrontFacets(
        opp_f_edge,
        opp_f_facets[0], opp_f_facets[1]);

    unique_ptr<FrontFacet3>* main_f_facets[3];
    unique_ptr<Vertex3>* main_vert;
    findAdjacentFrontFacets(
        frontEdge,
        main_f_facets[0], main_f_facets[1]);
    if ((*opp_f_facets[0])->facet->contains(**(*frontEdge)->edge->vertexes[0]))
    {
        main_f_facets[2] = opp_f_facets[0];
        main_vert = (*frontEdge)->edge->vertexes[0];
    }
    else if ((*opp_f_facets[0])->facet->contains(**(*frontEdge)->edge->vertexes[1]))
    {
        main_f_facets[2] = opp_f_facets[0];
        main_vert = (*frontEdge)->edge->vertexes[1];
    }
    else if ((*opp_f_facets[1])->facet->contains(**(*frontEdge)->edge->vertexes[0]))
    {
        main_f_facets[2] = opp_f_facets[1];
        main_vert = (*frontEdge)->edge->vertexes[0];
    }
    else
    {
        main_f_facets[2] = opp_f_facets[1];
        main_vert = (*frontEdge)->edge->vertexes[1];
    }

    unique_ptr<FrontEdge3>* new_f_facet_f_edge;
    unique_ptr<Edge3>* new_simp_edges[3];
    for (int i = 0; i < 3; i++)
    {
        new_simp_edges[i] = (*main_f_facets[i])->facet->findEdgeNotContainingVertex(**main_vert);
        new_f_facet_f_edge = findFrontEdge(new_simp_edges[i]);
        (*new_f_facet_f_edge)->needProcessing = true;
    }

    auto new_f_facet = addFrontFacet3(
        new_simp_edges[0],
        new_simp_edges[1],
        new_simp_edges[2]);
    Vec3 center = (*new_f_facet)->computeCenter();
    Vec3 opposite = (*main_vert)->getPosition();
    Vec3 third_pos = (*(*new_f_facet)->facet->findVertexNotIncludedInEdge(*(*oppositeEdge)->edge))->getPosition();
    Vec3 normal = Vec3::crossProduct(
        (*(*oppositeEdge)->edge->vertexes[0])->getPosition() - third_pos,
        (*(*oppositeEdge)->edge->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    (*new_f_facet)->setNormal(normal);

    innerSimps.push_back((new Simplex3(
        **(*frontEdge)->edge->vertexes[0],
        **(*frontEdge)->edge->vertexes[1],
        **(*opp_f_edge)->edge->vertexes[0],
        **(*opp_f_edge)->edge->vertexes[1]))->getPtrToUPtr());

    unique_ptr<FrontEdge3>* f_edge_ = nullptr;
    for (auto &f_facet : main_f_facets)
        for (auto &edge : (*f_facet)->facet->edges)
            if ((*edge)->contains(**main_vert) &&
                (f_edge_ = findFrontEdge(edge)))
                delete f_edge_->release();

    for (auto &f_facet : main_f_facets)
        delete f_facet->release();
}

void Crystallite3::exhaustWithoutNewVertexOppositeEdgeDontExists(unique_ptr<FrontEdge3>* frontEdge)
{
    unique_ptr<FrontEdge3>* opp_f_edge = nullptr;

    unique_ptr<FrontFacet3>* adj_f_facets[2];
    findAdjacentFrontFacets(
        frontEdge,
        adj_f_facets[0], adj_f_facets[1]);

    unique_ptr<Vertex3>* opp_verts[2];
    opp_verts[0] = (*adj_f_facets[0])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);
    opp_verts[1] = (*adj_f_facets[1])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);
    
    opp_f_edge = addFrontEdge3(opp_verts[0], opp_verts[1]);

    unique_ptr<Edge3>* new_simp_edges[3];
    new_simp_edges[2] = (*opp_f_edge)->edge->getPtrToUPtr();

    for (auto &edge : (*adj_f_facets[0])->facet->edges)
    {
        if ((*edge)->contains(**(*frontEdge)->edge->vertexes[0]) &&
            ((*edge)->contains(**(*opp_f_edge)->edge->vertexes[0]) ||
            (*edge)->contains(**(*opp_f_edge)->edge->vertexes[1])))
        {
            new_simp_edges[0] = edge;
            (*findFrontEdge(edge))->needProcessing = true;
            break;
        }
    }
    for (auto &edge : (*adj_f_facets[1])->facet->edges)
    {
        if ((*edge)->contains(**(*frontEdge)->edge->vertexes[0]) &&
            ((*edge)->contains(**(*opp_f_edge)->edge->vertexes[0]) ||
            (*edge)->contains(**(*opp_f_edge)->edge->vertexes[1])))
        {
            new_simp_edges[1] = edge;
            //if (!new_simp_edges[1] ||
            //    !*new_simp_edges[1])
            //{
            //    std::cout << (*new_simp_edges[1])->magnitude();
            //}
            (*findFrontEdge(edge))->needProcessing = true;
            break;
        }
    }
    auto new_f_facet = addFrontFacet3(
        new_simp_edges[0],
        new_simp_edges[1],
        new_simp_edges[2]);
    Vec3 center = (*new_f_facet)->computeCenter();
    Vec3 opposite = (*new_simp_edges[0])->contains(**(*frontEdge)->edge->vertexes[0]) ?
        (*(*frontEdge)->edge->vertexes[1])->getPosition() :
        (*(*frontEdge)->edge->vertexes[0])->getPosition();
    Vec3 third_pos = (*(*new_f_facet)->facet->findVertexNotIncludedInEdge(*(*opp_f_edge)->edge))->getPosition();
    Vec3 normal = Vec3::crossProduct(
        (*(*opp_f_edge)->edge->vertexes[0])->getPosition() - third_pos,
        (*(*opp_f_edge)->edge->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    (*new_f_facet)->setNormal(normal);

    for (auto &edge : (*adj_f_facets[0])->facet->edges)
        if ((*edge)->contains(**(*frontEdge)->edge->vertexes[1]) &&
            ((*edge)->contains(**(*opp_f_edge)->edge->vertexes[0]) ||
            (*edge)->contains(**(*opp_f_edge)->edge->vertexes[1])))
        {
            new_simp_edges[0] = edge;
            (*findFrontEdge(edge))->needProcessing = true;
            break;
        }
    for (auto &edge : (*adj_f_facets[1])->facet->edges)
        if ((*edge)->contains(**(*frontEdge)->edge->vertexes[1]) &&
            ((*edge)->contains(**(*opp_f_edge)->edge->vertexes[0]) ||
            (*edge)->contains(**(*opp_f_edge)->edge->vertexes[1])))
        {
            new_simp_edges[1] = edge;
            (*findFrontEdge(edge))->needProcessing = true;
            break;
        }
    new_f_facet = addFrontFacet3(
        new_simp_edges[0],
        new_simp_edges[1],
        new_simp_edges[2]);
    center = (*new_f_facet)->computeCenter();
    opposite = (*new_simp_edges[0])->contains(**(*frontEdge)->edge->vertexes[0]) ?
        (*(*frontEdge)->edge->vertexes[1])->getPosition() :
        (*(*frontEdge)->edge->vertexes[0])->getPosition();
    third_pos = (*(*new_f_facet)->facet->findVertexNotIncludedInEdge(*(*opp_f_edge)->edge))->getPosition();
    normal = Vec3::crossProduct(
        (*(*opp_f_edge)->edge->vertexes[0])->getPosition() - third_pos,
        (*(*opp_f_edge)->edge->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    (*new_f_facet)->setNormal(normal);

    innerSimps.push_back((new Simplex3(
        **(*frontEdge)->edge->vertexes[0],
        **(*frontEdge)->edge->vertexes[1],
        **opp_verts[0],
        **opp_verts[1]))->getPtrToUPtr());

    delete adj_f_facets[0]->release();
    delete adj_f_facets[1]->release();
    delete frontEdge->release();
}

void Crystallite3::exhaustWithoutNewVertex(unique_ptr<FrontEdge3>* frontEdge, const bool oppositeEdgeExistence, unique_ptr<FrontEdge3>* oppositeEdge)
{
    unique_ptr<FrontEdge3>* opp_f_edge = nullptr;
    if (oppositeEdgeExistence && oppositeEdge)
        opp_f_edge = oppositeEdge;
    else if (oppositeEdgeExistence)
        opp_f_edge = findOppositeFrontEdge(frontEdge);

    if ((oppositeEdge && oppositeEdgeExistence) ||
        opp_f_edge)
    {
        exhaustWithoutNewVertexOppositeEdgeExists(frontEdge, opp_f_edge);
        //std::cout << "__00__";
    }
    else
    {
        exhaustWithoutNewVertexOppositeEdgeDontExists(frontEdge);
        //std::cout << "__01__";
    }
}

const bool Crystallite3::tryComputeNewVertexPositionType3(unique_ptr<FrontFacet3>* frontFacet, Vec3& out_pos)
{
    unique_ptr<Edge3>* main_edges[2];
    main_edges[0] = (*frontFacet)->facet->edges[0];
    main_edges[1] = (*frontFacet)->facet->edges[1];
    unique_ptr<FrontEdge3>* main_f_edges[2];
    main_f_edges[0] = findFrontEdge(main_edges[0]);
    main_f_edges[1] = findFrontEdge(main_edges[1]);
    auto main_vert = (*main_edges[0])->vertexes[0] == (*main_edges[1])->vertexes[0] ?
        (*main_edges[0])->vertexes[0] :
        ((*main_edges[0])->vertexes[0] == (*main_edges[1])->vertexes[1] ?
        (*main_edges[0])->vertexes[0] :
            (*main_edges[0])->vertexes[1]);
    auto third_edge = (*frontFacet)->facet->findEdgeNotContainingVertex(**main_vert);
    auto third_f_edge = findFrontEdge(third_edge);

    double av_magn = 0.3333333333333333 * (
        (*(*frontFacet)->facet->edges[0])->magnitude()
        + (*(*frontFacet)->facet->edges[1])->magnitude()
        + (*(*frontFacet)->facet->edges[2])->magnitude());

    auto v0 = (*main_edges[0])->vertexes[0] == main_vert ? (*main_edges[0])->vertexes[1] : (*main_edges[0])->vertexes[0];
    auto v1 = (*main_edges[1])->vertexes[0] == main_vert ? (*main_edges[1])->vertexes[1] : (*main_edges[1])->vertexes[0];
    auto v2 = main_vert;

    unique_ptr<FrontFacet3>* adj_f_facets[2];
    findAdjacentFrontFacets(main_f_edges[0], adj_f_facets[0], adj_f_facets[1]);
    auto fn0 = adj_f_facets[0] == frontFacet ? adj_f_facets[1] : adj_f_facets[0];
    findAdjacentFrontFacets(main_f_edges[1], adj_f_facets[0], adj_f_facets[1]);
    auto fn1 = adj_f_facets[0] == frontFacet ? adj_f_facets[1] : adj_f_facets[0];
    findAdjacentFrontFacets(third_f_edge, adj_f_facets[0], adj_f_facets[1]);
    auto fn2 = adj_f_facets[0] == frontFacet ? adj_f_facets[1] : adj_f_facets[0];

    const Vec3& v0_pos = (*v0)->getPosition();
    const Vec3& v1_pos = (*v1)->getPosition();
    const Vec3& v2_pos = (*v2)->getPosition();

    Vec3 n_m = (*frontFacet)->getNormal();
    Vec3 n_n0 = (*fn0)->getNormal();
    Vec3 n_n1 = (*fn1)->getNormal();
    Vec3 n_n2 = (*fn2)->getNormal();

    Vec3 e_mn0 = (n_m + n_n0).normalize();
    Vec3 e_mn1 = (n_m + n_n1).normalize();
    Vec3 e_mn2 = n_m + n_n2;

    Vec3 np_mn0 = Vec3::crossProduct(v2_pos - v0_pos, e_mn0);
    Vec3 np_mn1 = Vec3::crossProduct(v2_pos - v1_pos, e_mn1);

    Vec3 e = Vec3::crossProduct(np_mn0, np_mn1);

    Vec3 new_pos = spatialalgs::lineIntersectPlane(v2_pos, e, v0_pos, v1_pos, v0_pos + e_mn2);

    Vec3 v0_to_np = new_pos - v0_pos;
    Vec3 v1_to_np = new_pos - v1_pos;
    Vec3 v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + 1e-3 * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + 1e-3 * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + 1e-3 * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
        !vertexInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

    out_pos = new_pos;
    return true;
}

const bool Crystallite3::tryComputeNewVertexPositionType2(const int smallAngleIndex0, const double angleCos0, const int smallAngleIndex1, const double angleCos1, unique_ptr<FrontFacet3>* frontFacet, Vec3& out_pos)
{
    unique_ptr<Edge3>* main_edges[2];
    main_edges[0] = (*frontFacet)->facet->edges[smallAngleIndex0];
    main_edges[1] = (*frontFacet)->facet->edges[smallAngleIndex1];
    unique_ptr<FrontEdge3>* main_f_edges[2];
    main_f_edges[0] = findFrontEdge(main_edges[0]);
    main_f_edges[1] = findFrontEdge(main_edges[1]);
    auto main_vert = (*main_edges[0])->vertexes[0] == (*main_edges[1])->vertexes[0] ?
        (*main_edges[0])->vertexes[0] :
        ((*main_edges[0])->vertexes[0] == (*main_edges[1])->vertexes[1] ?
            (*main_edges[0])->vertexes[0] :
            (*main_edges[0])->vertexes[1]);

    double av_magn = 0.3333333333333333 * (
        (*(*frontFacet)->facet->edges[0])->magnitude()
        + (*(*frontFacet)->facet->edges[1])->magnitude()
        + (*(*frontFacet)->facet->edges[2])->magnitude());

    auto v0 = (*main_edges[0])->vertexes[0] == main_vert ? (*main_edges[0])->vertexes[1] : (*main_edges[0])->vertexes[0];
    auto v1 = (*main_edges[1])->vertexes[0] == main_vert ? (*main_edges[1])->vertexes[1] : (*main_edges[1])->vertexes[0];
    auto v2 = main_vert;
    unique_ptr<FrontFacet3>* adj_f_facets[2];
    findAdjacentFrontFacets(main_f_edges[0], adj_f_facets[0], adj_f_facets[1]);
    auto fn0 = adj_f_facets[0] == frontFacet ? adj_f_facets[1] : adj_f_facets[0];
    findAdjacentFrontFacets(main_f_edges[1], adj_f_facets[0], adj_f_facets[1]);
    auto fn1 = adj_f_facets[0] == frontFacet ? adj_f_facets[1] : adj_f_facets[0];

    const Vec3& v0_pos = (*v0)->getPosition();
    const Vec3& v1_pos = (*v1)->getPosition();
    const Vec3& v2_pos = (*v2)->getPosition();

    Vec3 n_m = (*frontFacet)->getNormal();
    Vec3 n_n0 = (*fn0)->getNormal();
    Vec3 n_n1 = (*fn1)->getNormal();
    
    Vec3 e_mn0 = (n_m + n_n0).normalize();
    Vec3 e_mn1 = (n_m + n_n1).normalize();

    Vec3 np_mn0 = Vec3::crossProduct(v2_pos - v0_pos, e_mn0);
    Vec3 np_mn1 = Vec3::crossProduct(v2_pos - v1_pos, e_mn1);

    Vec3 e = Vec3::crossProduct(np_mn0, np_mn1).normalize();
    if (Vec3::dotProduct(e, n_m)) e *= -1.0;

    Vec3 new_pos = angleCos0 > COS_DEG_80 && angleCos1 > COS_DEG_80 ?
        v2_pos + 0.5 * av_magn * e :
        (angleCos0 > COS_DEG_80 || angleCos1 > COS_DEG_80 ?
            v2_pos + 0.7 * av_magn * e :
            v2_pos + 0.9 * av_magn * e);

    Vec3 v0_to_np = new_pos - v0_pos;
    Vec3 v1_to_np = new_pos - v1_pos;
    Vec3 v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + 1e-3 * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + 1e-3 * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + 1e-3 * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
        !vertexInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

    out_pos = new_pos;
    return true;
}

const bool Crystallite3::tryComputeNewVertexPositionType1(const int smallAngleIndex, const double angleCos, unique_ptr<FrontFacet3>* frontFacet, Vec3& out_pos)
{
    auto main_edge = (*frontFacet)->facet->edges[smallAngleIndex];
    auto main_f_edge = findFrontEdge(main_edge);
    Vec3 main_edge_poses[2];
    main_edge_poses[0] = (*(*main_edge)->vertexes[0])->getPosition();
    main_edge_poses[1] = (*(*main_edge)->vertexes[1])->getPosition();

    double av_magn = 0.3333333333333333 * (
        (*(*frontFacet)->facet->edges[0])->magnitude()
        + (*(*frontFacet)->facet->edges[1])->magnitude()
        + (*(*frontFacet)->facet->edges[2])->magnitude());

    auto v0 = (*main_edge)->vertexes[0];
    auto v1 = (*main_edge)->vertexes[1];
    auto v2 = (*frontFacet)->facet->findVertexNotIncludedInEdge(**main_edge);
    unique_ptr<FrontFacet3>* adj_f_facets[2];
    findAdjacentFrontFacets(main_f_edge, adj_f_facets[0], adj_f_facets[1]);
    auto fn = adj_f_facets[0] == frontFacet ? adj_f_facets[1] : adj_f_facets[0];
    auto vn = (*fn)->facet->findVertexNotIncludedInEdge(**main_edge);

    const Vec3& v0_pos = (*v0)->getPosition();
    const Vec3& v1_pos = (*v1)->getPosition();
    const Vec3& v2_pos = (*v2)->getPosition();
    const Vec3& vn_pos = (*vn)->getPosition();

    Vec3 v2pr = spatialalgs::project(v2_pos, v0_pos, v1_pos);
    Vec3 vnpr = spatialalgs::project(vn_pos, v0_pos, v1_pos);

    Vec3 c = 0.25 * (v0_pos + v1_pos + v2pr + vnpr);
    Vec3 e = ((*frontFacet)->getNormal() + (*fn)->getNormal()).normalize();
    Vec3 new_pos = angleCos > COS_DEG_80 ?
        c + 0.6 * ONE_PLUS_SQRT2_DIV_SQRT3 * av_magn * e :
        c + 0.9 * ONE_PLUS_SQRT2_DIV_SQRT3 * av_magn * e;

    Vec3 v0_to_np = new_pos - v0_pos;
    Vec3 v1_to_np = new_pos - v1_pos;
    Vec3 v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + 1e-3 * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + 1e-3 * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + 1e-3 * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
        !vertexInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

    out_pos = new_pos;
    return true;
}

const bool Crystallite3::tryComputeNewVertexPositionType0(unique_ptr<FrontFacet3>* frontFacet, Vec3& out_pos)
{
    double av_magn = 0.3333333333333333 * (
          (*(*frontFacet)->facet->edges[0])->magnitude()
        + (*(*frontFacet)->facet->edges[1])->magnitude()
        + (*(*frontFacet)->facet->edges[2])->magnitude());
    Vec3 new_pos = (*frontFacet)->computeCenter() + SQRT_2_DIV_3 * av_magn * (*frontFacet)->getNormal();

    auto v0 = (*(*frontFacet)->facet->edges[0])->vertexes[0];
    auto v1 = (*(*frontFacet)->facet->edges[0])->vertexes[1];
    auto v2 = (*frontFacet)->facet->findVertexNotIncludedInEdge(**(*frontFacet)->facet->edges[0]);
    const Vec3& v0_pos = (*v0)->getPosition();
    const Vec3& v1_pos = (*v1)->getPosition();
    const Vec3& v2_pos = (*v2)->getPosition();
    Vec3 v0_to_np = new_pos - v0_pos;
    Vec3 v1_to_np = new_pos - v1_pos;
    Vec3 v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + 1e-3 * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + 1e-3 * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + 1e-3 * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
        !vertexInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

    out_pos = new_pos;
    return true;
}

const bool Crystallite3::tryComputeNewVertexPosition(unique_ptr<FrontFacet3>* frontFacet, Vec3& out_pos)
{
    double angs_coses[3]
    { 
        (*findFrontEdge((*frontFacet)->facet->edges[0]))->getAngleExCos(this),
        (*findFrontEdge((*frontFacet)->facet->edges[1]))->getAngleExCos(this),
        (*findFrontEdge((*frontFacet)->facet->edges[2]))->getAngleExCos(this)
    };
    int indexes[3];
    int small_angs_num = 0;
    if (angs_coses[0] > COS_DEG_140) indexes[small_angs_num++] = 0;
    if (angs_coses[1] > COS_DEG_140) indexes[small_angs_num++] = 1;
    if (angs_coses[2] > COS_DEG_140) indexes[small_angs_num++] = 2;

    switch (small_angs_num)
    {
    case 0: return tryComputeNewVertexPositionType0(frontFacet, out_pos);        
    case 1: return tryComputeNewVertexPositionType1(indexes[0], angs_coses[indexes[0]], frontFacet, out_pos);
    case 2: return tryComputeNewVertexPositionType2(indexes[0], angs_coses[indexes[0]], indexes[1], angs_coses[indexes[1]], frontFacet, out_pos);
    case 3: return tryComputeNewVertexPositionType3(frontFacet, out_pos);
    }

    return true;
}

unique_ptr<FrontFacet3>* Crystallite3::chooseFrontFacetForExhaustionWithNewVertex(unique_ptr<FrontEdge3>* frontEdge)
{
    unique_ptr<FrontFacet3>* adj_f_facets[2];
    findAdjacentFrontFacets(
        frontEdge,
        adj_f_facets[0], adj_f_facets[1]);

    return (*adj_f_facets[0])->computeQuality() < (*adj_f_facets[1])->computeQuality() ?
        adj_f_facets[0] :
        adj_f_facets[1];
}

void Crystallite3::exhaustWithNewVertex(unique_ptr<FrontFacet3>* frontFacet, Vec3 vertPos)
{
    unique_ptr<Vertex3>* new_vert = (new Vertex3(vertPos))->getPtrToUPtr();
    innerVerts.push_back(new_vert);

    unique_ptr<Edge3>* new_simp_edges[6];
    new_simp_edges[0] = (*(*frontFacet)->facet->edges[0])->getPtrToUPtr();
    auto far_vert = (*frontFacet)->facet->findVertexNotIncludedInEdge(**new_simp_edges[0]);
    new_simp_edges[1] = (*frontFacet)->facet->findEdge(**(*new_simp_edges[0])->vertexes[1], **far_vert);
    new_simp_edges[2] = (*frontFacet)->facet->findEdge(**(*new_simp_edges[0])->vertexes[0], **far_vert);
    new_simp_edges[3] = (*addFrontEdge3((*new_simp_edges[0])->vertexes[0], new_vert))->edge->getPtrToUPtr();
    new_simp_edges[4] = (*addFrontEdge3((*new_simp_edges[0])->vertexes[1], new_vert))->edge->getPtrToUPtr();
    new_simp_edges[5] = (*addFrontEdge3(far_vert, new_vert))->edge->getPtrToUPtr();

    (*findFrontEdge(new_simp_edges[0]))->needProcessing = true;
    (*findFrontEdge(new_simp_edges[1]))->needProcessing = true;
    (*findFrontEdge(new_simp_edges[2]))->needProcessing = true;
    
    auto new_f_facet = addFrontFacet3(
        new_simp_edges[0],
        new_simp_edges[3],
        new_simp_edges[4]);
    (*new_f_facet)->setNormal(computeNormalInSimplex3(new_f_facet, (*far_vert)->getPosition()));
    //std::cout << " | " << (*new_f_facet)->getNormal().coors[0] 
    //    << " " << (*new_f_facet)->getNormal().coors[1]
    //    << " " << (*new_f_facet)->getNormal().coors[2] << " | ";

    new_f_facet = addFrontFacet3(
        new_simp_edges[2],
        new_simp_edges[3],
        new_simp_edges[5]);
    (*new_f_facet)->setNormal(computeNormalInSimplex3(new_f_facet, (*(*new_simp_edges[0])->vertexes[1])->getPosition()));
    //std::cout << " | " << (*new_f_facet)->getNormal().coors[0]
    //    << " " << (*new_f_facet)->getNormal().coors[1]
    //    << " " << (*new_f_facet)->getNormal().coors[2] << " | ";

    new_f_facet = addFrontFacet3(
        new_simp_edges[1],
        new_simp_edges[4],
        new_simp_edges[5]);
    (*new_f_facet)->setNormal(computeNormalInSimplex3(new_f_facet, (*(*new_simp_edges[0])->vertexes[0])->getPosition()));
    //std::cout << " | " << (*new_f_facet)->getNormal().coors[0]
    //    << " " << (*new_f_facet)->getNormal().coors[1]
    //    << " " << (*new_f_facet)->getNormal().coors[2] << " | ";

    innerSimps.push_back((new Simplex3(
        **(*new_simp_edges[0])->vertexes[0],
        **(*new_simp_edges[0])->vertexes[1],
        **(*new_simp_edges[5])->vertexes[0],
        **(*new_simp_edges[5])->vertexes[1]))->getPtrToUPtr());

    delete frontFacet->release();
}

const bool Crystallite3::NewVertexPosition_OLD(unique_ptr<FrontEdge3>* frontEdge, Vec3& out_pos)
{
    Vec3 f_edge_verts_poses[2];
    f_edge_verts_poses[0] = (*(*frontEdge)->edge->vertexes[0])->getPosition();
    f_edge_verts_poses[1] = (*(*frontEdge)->edge->vertexes[1])->getPosition();

    Vec3 orig = 0.5 * (f_edge_verts_poses[0] + f_edge_verts_poses[1]);

    unique_ptr<FrontFacet3>* f_facets[2];
    findAdjacentFrontFacets(frontEdge, f_facets[0], f_facets[1]);

    unique_ptr<Vertex3>* opp_verts[2];
    opp_verts[0] = (*f_facets[0])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);
    opp_verts[1] = (*f_facets[1])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);

    Vec3 opp_verts_poses[2];
    opp_verts_poses[0] = (*opp_verts[0])->getPosition();
    opp_verts_poses[1] = (*opp_verts[1])->getPosition();

    Vec3 facets_inner_vecs[2];
    facets_inner_vecs[0] = opp_verts_poses[0] - orig;
    facets_inner_vecs[1] = opp_verts_poses[1] - orig;

    Vec3 dir = (facets_inner_vecs[0] + facets_inner_vecs[1]).normalize() * _preferredLength * 0.866; // sind(60) == 0.8660...
    Vec3 orig_plus_ex_dir = orig + dir * 1.6;
    while (
        segmentFrontIntersectionCheck(f_edge_verts_poses[0] + 1e-3 * (orig_plus_ex_dir - f_edge_verts_poses[0]), orig_plus_ex_dir) ||
        segmentFrontIntersectionCheck(f_edge_verts_poses[1] + 1e-3 * (orig_plus_ex_dir - f_edge_verts_poses[1]), orig_plus_ex_dir) ||
        segmentFrontIntersectionCheck(opp_verts_poses[0] + 1e-3 * (orig_plus_ex_dir - opp_verts_poses[0]), orig_plus_ex_dir) ||
        segmentFrontIntersectionCheck(opp_verts_poses[1] + 1e-3 * (orig_plus_ex_dir - opp_verts_poses[1]), orig_plus_ex_dir) ||
        !vertexInsideFrontCheck(orig_plus_ex_dir) ||
        facetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[0], orig_plus_ex_dir) ||
        facetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], orig_plus_ex_dir) ||
        facetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[1], orig_plus_ex_dir) ||
        facetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[1], orig_plus_ex_dir)
        //!vertexInsideFrontCheck(orig + dir) ||
        /*lineSegmentGlobalIntersectionCheck(f_edge_verts_poses[0] + 1e-3 * (orig + dir - f_edge_verts_poses[0]), orig + dir) ||
        lineSegmentGlobalIntersectionCheck(f_edge_verts_poses[1] + 1e-3 * (orig + dir - f_edge_verts_poses[1]), orig + dir) ||
        lineSegmentGlobalIntersectionCheck(f_facets_far_verts_poses[0] + 1e-3 * (orig + dir - f_facets_far_verts_poses[0]), orig + dir) ||
        segmentFrontIntersectionCheck(f_facets_far_verts_poses[1] + 1e-3 * (orig + dir - f_facets_far_verts_poses[1]), orig + dir)*/)
    {
        dir *= 0.75;
        orig_plus_ex_dir = orig + dir * 1.6;

        if (dir.magnitude() < _preferredLength * 0.2)
            return false;
    }

    /*if (dir.magnitude() < _preferredLength * 0.03)
        return false;*/

    /*if (facetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[0], orig + dir) ||
        facetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], orig + dir) ||
        facetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[1], orig + dir) ||
        facetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], orig + dir))
        return false;*/

    /*for (auto &vert : innerVerts)
    {
        if (!*vert ||
            vert == (*frontEdge)->edge->vertexes[0] ||
            vert == (*frontEdge)->edge->vertexes[1] ||
            vert == opp_verts[0] ||
            vert == opp_verts[1])
            continue;

        if (insideSimplex3Check(
                (*(*frontEdge)->edge->vertexes[0])->getPosition(),
                (*(*frontEdge)->edge->vertexes[1])->getPosition(),
                (*opp_verts[0])->getPosition(),
                orig_plus_ex_dir,
                (*vert)->getPosition()) ||
            insideSimplex3Check(
                (*(*frontEdge)->edge->vertexes[0])->getPosition(),
                (*(*frontEdge)->edge->vertexes[1])->getPosition(),
                (*opp_verts[1])->getPosition(),
                orig_plus_ex_dir,
                (*vert)->getPosition()))
            return false;
    }*/

    out_pos = orig + dir;
    return true;
}

void Crystallite3::ExhaustWithNewVertex_OLD(unique_ptr<FrontEdge3>* frontEdge, Vec3 vertPos)
{
    unique_ptr<FrontFacet3>* adj_f_facets[2];
    findAdjacentFrontFacets(frontEdge, adj_f_facets[0], adj_f_facets[1]);

    unique_ptr<Vertex3>* opp_verts[2];
    opp_verts[0] = (*adj_f_facets[0])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);
    opp_verts[1] = (*adj_f_facets[1])->facet->findVertexNotIncludedInEdge(*(*frontEdge)->edge);

    unique_ptr<Vertex3>* new_vert = (new Vertex3(vertPos))->getPtrToUPtr();
    innerVerts.push_back(new_vert);

    unique_ptr<Edge3>* new_simp_edges[9];
    new_simp_edges[0] = (*frontEdge)->edge->getPtrToUPtr();
    new_simp_edges[1] = (*addFrontEdge3((*frontEdge)->edge->vertexes[0], new_vert))->edge->getPtrToUPtr();
    new_simp_edges[2] = (*addFrontEdge3((*frontEdge)->edge->vertexes[1], new_vert))->edge->getPtrToUPtr();
    new_simp_edges[3] = (*adj_f_facets[0])->facet->findEdge(**(*frontEdge)->edge->vertexes[0], **opp_verts[0]);
    new_simp_edges[4] = (*addFrontEdge3(opp_verts[0], new_vert))->edge->getPtrToUPtr();
    new_simp_edges[5] = (*adj_f_facets[0])->facet->findEdge(**(*frontEdge)->edge->vertexes[1], **opp_verts[0]);
    new_simp_edges[6] = (*adj_f_facets[1])->facet->findEdge(**(*frontEdge)->edge->vertexes[0], **opp_verts[1]);
    new_simp_edges[7] = (*addFrontEdge3(opp_verts[1], new_vert))->edge->getPtrToUPtr();
    new_simp_edges[8] = (*adj_f_facets[1])->facet->findEdge(**(*frontEdge)->edge->vertexes[1], **opp_verts[1]);

    (*findFrontEdge(new_simp_edges[3]))->needProcessing = true;
    (*findFrontEdge(new_simp_edges[5]))->needProcessing = true;
    (*findFrontEdge(new_simp_edges[6]))->needProcessing = true;
    (*findFrontEdge(new_simp_edges[8]))->needProcessing = true;

    innerFacets.push_back((new Facet3(
        **new_simp_edges[0], 
        **new_simp_edges[1], 
        **new_simp_edges[2]))->getPtrToUPtr());

    unique_ptr<FrontFacet3>* new_f_facet;
    new_f_facet = addFrontFacet3(
        new_simp_edges[1],
        new_simp_edges[3],
        new_simp_edges[4]);
    Vec3 center = (*new_f_facet)->computeCenter();
    Vec3 opposite = (*new_f_facet)->facet->contains(**(*new_simp_edges[0])->vertexes[0]) ?
        (*(*new_simp_edges[0])->vertexes[1])->getPosition() :
        (*(*new_simp_edges[0])->vertexes[0])->getPosition();
    Vec3 third_pos = (*(*new_f_facet)->facet->findVertexNotIncludedInEdge(**new_simp_edges[4]))->getPosition();
    Vec3 normal = Vec3::crossProduct(
        (*(*new_simp_edges[4])->vertexes[0])->getPosition() - third_pos,
        (*(*new_simp_edges[4])->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    (*new_f_facet)->setNormal(normal);

    new_f_facet = addFrontFacet3(
        new_simp_edges[2],
        new_simp_edges[4],
        new_simp_edges[5]);
    center = (*new_f_facet)->computeCenter();
    opposite = (*new_f_facet)->facet->contains(**(*new_simp_edges[0])->vertexes[0]) ?
        (*(*new_simp_edges[0])->vertexes[1])->getPosition() :
        (*(*new_simp_edges[0])->vertexes[0])->getPosition();
    third_pos = (*(*new_f_facet)->facet->findVertexNotIncludedInEdge(**new_simp_edges[4]))->getPosition();
    normal = Vec3::crossProduct(
        (*(*new_simp_edges[4])->vertexes[0])->getPosition() - third_pos,
        (*(*new_simp_edges[4])->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    (*new_f_facet)->setNormal(normal);

    new_f_facet = addFrontFacet3(
        new_simp_edges[1],
        new_simp_edges[6],
        new_simp_edges[7]);
    center = (*new_f_facet)->computeCenter();
    opposite = (*new_f_facet)->facet->contains(**(*new_simp_edges[0])->vertexes[0]) ?
        (*(*new_simp_edges[0])->vertexes[1])->getPosition() :
        (*(*new_simp_edges[0])->vertexes[0])->getPosition();
    third_pos = (*(*new_f_facet)->facet->findVertexNotIncludedInEdge(**new_simp_edges[7]))->getPosition();
    normal = Vec3::crossProduct(
        (*(*new_simp_edges[7])->vertexes[0])->getPosition() - third_pos,
        (*(*new_simp_edges[7])->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    (*new_f_facet)->setNormal(normal);

    new_f_facet = addFrontFacet3(
        new_simp_edges[2],
        new_simp_edges[8],
        new_simp_edges[7]);
    center = (*new_f_facet)->computeCenter();
    opposite = (*new_f_facet)->facet->contains(**(*new_simp_edges[0])->vertexes[0]) ?
        (*(*new_simp_edges[0])->vertexes[1])->getPosition() :
        (*(*new_simp_edges[0])->vertexes[0])->getPosition();
    third_pos = (*(*new_f_facet)->facet->findVertexNotIncludedInEdge(**new_simp_edges[7]))->getPosition();
    normal = Vec3::crossProduct(
        (*(*new_simp_edges[7])->vertexes[0])->getPosition() - third_pos,
        (*(*new_simp_edges[7])->vertexes[1])->getPosition() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    (*new_f_facet)->setNormal(normal);

    innerSimps.push_back((new Simplex3(
        **(*new_simp_edges[0])->vertexes[0],
        **(*new_simp_edges[0])->vertexes[1],
        **(*new_simp_edges[4])->vertexes[0],
        **(*new_simp_edges[4])->vertexes[1]))->getPtrToUPtr());
    innerSimps.push_back((new Simplex3(
        **(*new_simp_edges[0])->vertexes[0],
        **(*new_simp_edges[0])->vertexes[1],
        **(*new_simp_edges[7])->vertexes[0],
        **(*new_simp_edges[7])->vertexes[1]))->getPtrToUPtr());

    delete adj_f_facets[0]->release();
    delete adj_f_facets[1]->release();
    delete frontEdge->release();
}

const bool Crystallite3::tryExhaustWithoutNewVertex(unique_ptr<FrontEdge3>* frontEdge, const bool oppositeEdgeExistence, unique_ptr<FrontEdge3>* oppositeEdge)
{
    if (parallelFacetsCheck(frontEdge) ||
        edgeIntersectionCheck(frontEdge) ||
        facetsIntersectionCheck(frontEdge) ||
        anyVertexInsidePotentialSimplex3Check(frontEdge) ||
        frontSplitCheck(frontEdge))
    {
        return false;
    }

    exhaustWithoutNewVertex(frontEdge, oppositeEdgeExistence, oppositeEdge);
    return true;
}

const bool Crystallite3::tryExhaustWithNewVertex(unique_ptr<FrontEdge3>* frontEdge)
{
    if (edgeIntersectionCheck(frontEdge) ||
        parallelFacetsCheck(frontEdge))
    {
        return false;
    }

    auto exhaust_f_facet = chooseFrontFacetForExhaustionWithNewVertex(frontEdge);
    Vec3 new_vert_pos;
    if (!tryComputeNewVertexPosition(exhaust_f_facet, new_vert_pos))
        return false;

    exhaustWithNewVertex(exhaust_f_facet, new_vert_pos);
    return true;
}

template<class T>
size_t PtrsToNullptrNumber(vector<unique_ptr<T>*> &vec)
{
    return std::count_if(vec.begin(), vec.end(), [](unique_ptr<T>* &ptr) { return !*ptr ? true : false; });
}

const bool Crystallite3::ProcessVerySmallAngles(Polycrystal3* polycr)
{
    if (frontContainsOfOnly1Simplex3OrEmpty())
    {
        // Do something...

        for (auto &f_facet : frontFacets)
        {
            if (*f_facet)
                delete f_facet->release();
        }
        for (auto &f_edge : frontEdges)
        {
            if (*f_edge)
                delete f_edge->release();
        }
        return true;
    }

    for (size_t i = 0ull; i < frontEdges.size();)
    {
        if (!*frontEdges[i])
        {
            i++;
            continue;
        }

        if ((*frontEdges[i])->getAngleExCos(this) <= COS_DEG_80 ||
            parallelFacetsCheck(frontEdges[i]) ||
            edgeIntersectionCheck(frontEdges[i]) ||
            facetsIntersectionCheck(frontEdges[i]) ||
            anyVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
            frontSplitCheck(frontEdges[i]))
        {
            i++;
            continue;
        }

        /*if (!*frontEdges[i] || 
            (*frontEdges[i])->computeAngle(frontFacets) >= DEG_80_IN_RADIANS ||
            parallelFacetsCheck(frontEdges[i]) ||
            edgeIntersectionCheck(frontEdges[i]) ||
            facetsIntersectionCheck(frontEdges[i]) ||
            anyVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
            frontSplitCheck(frontEdges[i]))
        {
            i++;
            continue;
        }*/

        exhaustWithoutNewVertex(frontEdges[i]);
        i = 0ull;

        /*if (globalIntersectionCheck())
        std::cout << "Bad!";*/
        //polycr->outputData();

        if (frontContainsOfOnly1Simplex3OrEmpty())
        {
            // Do something...

            for (auto &f_facet : frontFacets)
            {
                if (*f_facet)
                    delete f_facet->release();
            }
            for (auto &f_edge : frontEdges)
            {
                if (*f_edge)
                    delete f_edge->release();
            }
            return true;
        }

        if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
            removePtrsToNullptr(frontEdges);
        if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
            removePtrsToNullptr(frontFacets);
    }

    return false;
}

const bool Crystallite3::ProcessSmallAngles(Polycrystal3* polycr)
{
    if (frontContainsOfOnly1Simplex3OrEmpty())
    {
        // Do something...

        for (auto &f_facet : frontFacets)
        {
            if (*f_facet)
                delete f_facet->release();
        }
        for (auto &f_edge : frontEdges)
        {
            if (*f_edge)
                delete f_edge->release();
        }
        return true;
    }

    for (size_t i = 0ull; i < frontEdges.size();)
    {
        if (!*frontEdges[i])
        {
            i++;
            continue;
        }

        if ((*frontEdges[i])->getAngleExCos(this) <= COS_DEG_100 ||
            parallelFacetsCheck(frontEdges[i]) ||
            edgeIntersectionCheck(frontEdges[i]) ||
            facetsIntersectionCheck(frontEdges[i]) ||
            anyVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
            frontSplitCheck(frontEdges[i]))
        {
            i++;
            continue;
        }

        /*if (!*frontEdges[i] ||
            (*frontEdges[i])->computeAngle(frontFacets) >= DEG_100_IN_RADIANS ||
            parallelFacetsCheck(frontEdges[i]) ||
            edgeIntersectionCheck(frontEdges[i]) ||
            facetsIntersectionCheck(frontEdges[i]) ||
            anyVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
            frontSplitCheck(frontEdges[i]))
        {
            i++;
            continue;
        }*/

        exhaustWithoutNewVertex(frontEdges[i]);
        i = 0ull;

        /*if (globalIntersectionCheck())
            std::cout << "Bad!";*/
        //polycr->outputData();

        if (frontContainsOfOnly1Simplex3OrEmpty())
        {
            // Do something...

            for (auto &f_facet : frontFacets)
            {
                if (*f_facet)
                    delete f_facet->release();
            }
            for (auto &f_edge : frontEdges)
            {
                if (*f_edge)
                    delete f_edge->release();
            }
            return true;
        }

        if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
            removePtrsToNullptr(frontEdges);
        if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
            removePtrsToNullptr(frontFacets);
    }

    return false;
}

const bool Crystallite3::ProcessMediumAngles(Polycrystal3* polycr)
{
    if (frontContainsOfOnly1Simplex3OrEmpty())
    {
        // Do something...

        for (auto &f_facet : frontFacets)
        {
            if (*f_facet)
                delete f_facet->release();
        }
        for (auto &f_edge : frontEdges)
        {
            if (*f_edge)
                delete f_edge->release();
        }
        return true;
    }

    for (size_t i = 0ull; i < frontEdges.size();)
    {
        if (!*frontEdges[i])
        {
            i++;
            continue;
        }

        Vec3 new_vert_pos;
        if ((*frontEdges[i])->getAngleExCos(this) <= COS_DEG_150 ||
            edgeIntersectionCheck(frontEdges[i]) ||
            parallelFacetsCheck(frontEdges[i]) ||
            !NewVertexPosition_OLD(frontEdges[i], new_vert_pos))
        {
            i++;
            continue;
        }

        /*Vec3 new_vert_pos;
        if (!*frontEdges[i] ||
            (*frontEdges[i])->computeAngle(frontFacets) >= DEG_150_IN_RADIANS ||
            edgeIntersectionCheck(frontEdges[i]) ||
            parallelFacetsCheck(frontEdges[i]) ||
            !NewVertexPosition_OLD(frontEdges[i], new_vert_pos))
        {
            i++;
            continue;
        }*/

        ExhaustWithNewVertex_OLD(frontEdges[i], new_vert_pos);
        if (ProcessVerySmallAngles(polycr))
            return true;
        if (ProcessSmallAngles(polycr))
            return true;
        //polycr->outputData();
        i = 0ull;

        if (frontContainsOfOnly1Simplex3OrEmpty())
        {
            // Do something...

            for (auto &f_facet : frontFacets)
            {
                if (*f_facet)
                    delete f_facet->release();
            }
            for (auto &f_edge : frontEdges)
            {
                if (*f_edge)
                    delete f_edge->release();
            }

            return true;
        }

        if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
            removePtrsToNullptr(frontEdges);
        if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
            removePtrsToNullptr(frontFacets);
    }

    return false;
}

void Crystallite3::ProcessLargeAngles(Polycrystal3* polycr)
{
    if (frontContainsOfOnly1Simplex3OrEmpty())
    {
        // Do something...

        for (auto &f_facet : frontFacets)
        {
            if (*f_facet)
                delete f_facet->release();
        }
        for (auto &f_edge : frontEdges)
        {
            if (*f_edge)
                delete f_edge->release();
        }
        return;
    }
    
    for (size_t i = 0ull; i < frontEdges.size();)
    {
        if (!*frontEdges[i])
        {
            i++;
            continue;
        }

        Vec3 new_vert_pos;
        if ((*frontEdges[i])->getAngleExCos(this) <= COS_DEG_180 ||
            edgeIntersectionCheck(frontEdges[i]) ||
            parallelFacetsCheck(frontEdges[i]) ||
            !NewVertexPosition_OLD(frontEdges[i], new_vert_pos))
        {
            i++;
            continue;
        }

        ExhaustWithNewVertex_OLD(frontEdges[i], new_vert_pos);
        /*if (globalIntersectionCheck())
            std::cout << "Bad!";*/
        //polycr->outputData();
        if (ProcessMediumAngles(polycr))
            return;
        /*if (globalIntersectionCheck())
            std::cout << "Bad!";*/
        //polycr->outputData();
        i = 0ull;

        if (frontContainsOfOnly1Simplex3OrEmpty())
        {
            // Do something...

            for (auto &f_facet : frontFacets)
            {
                if (*f_facet)
                    delete f_facet->release();
            }
            for (auto &f_edge : frontEdges)
            {
                if (*f_edge)
                    delete f_edge->release();
            }
            return;
        }

        if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
            removePtrsToNullptr(frontEdges);
        if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
            removePtrsToNullptr(frontFacets);
    }

    polycr->outputData();
    if (!frontContainsOfOnly1Simplex3OrEmpty())
        throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");
}

void Crystallite3::processAngles(Polycrystal3* polycr)
{
    if (frontContainsOfOnly1Simplex3OrEmpty())
    {
        throw std::logic_error("Error in function: Crystallite::processAngles");

        for (auto &f_facet : frontFacets)
        {
            if (*f_facet)
                delete f_facet->release();
        }
        for (auto &f_edge : frontEdges)
        {
            if (*f_edge)
                delete f_edge->release();
        }
        return;
    }

    double max_excos = 2.0; //int exhaustions = 0;
    for (unique_ptr<FrontEdge3>* curr_f_edge = currentFrontEdge(max_excos); ; curr_f_edge = currentFrontEdge(max_excos))
    {
        if (!curr_f_edge)
        {
            polycr->outputData();
            throw std::logic_error("Crystallite3::currentFrontEdge returned nullptr");
        }
        if (exhaustWithoutNewVertexPriorityPredicate(curr_f_edge))
        {
            if (!tryExhaustWithoutNewVertex(curr_f_edge))
            {
                max_excos = (*curr_f_edge)->getAngleExCos(this);
                continue;
            }
            //std::cout << "__0__";
        }
        else if (exhaustWithNewVertexPriorityPredicate(curr_f_edge))
        {
            if (!tryExhaustWithNewVertex(curr_f_edge))
            {
                max_excos = (*curr_f_edge)->getAngleExCos(this);
                continue;
            }
            //std::cout << "__1__";
        }
        else
        {
            unique_ptr<FrontFacet3>* exhaust_from_f_facet = nullptr;
            Vec3* new_vert_pos = nullptr;
            switch (exhaustTypeQualityPriorityCalculation(curr_f_edge, &exhaust_from_f_facet, &new_vert_pos))
            {
            case WITHOUT_NEW_VERTEX:
                exhaustWithoutNewVertex(curr_f_edge);
                //std::cout << "__2__";
                break;

            case WITH_NEW_VERTEX:
                if (new_vert_pos)
                {
                    exhaustWithNewVertex(exhaust_from_f_facet, *new_vert_pos);
                    delete new_vert_pos;
                    //delete exhaust_from_f_facet;
                    //std::cout << "__3__";
                }
                else
                {
                    if (!tryExhaustWithNewVertex(curr_f_edge))
                    {
                        max_excos = (*curr_f_edge)->getAngleExCos(this);
                        continue;
                    }
                    //std::cout << "__4__";
                }
                break;

            case DONT_EXHAUST:
                max_excos = (*curr_f_edge)->getAngleExCos(this);
                continue;
                break;
            }
        }
        max_excos = 2.0;
        //polycr->outputData();

        //std::cout << " exhs = " << ++exhaustions << "\n";
        //if (exhaustions >= 70)
        //    polycr->outputData();
        //debug();

        if (frontContainsOfOnly1Simplex3OrEmpty())
        {
            unique_ptr<FrontFacet3>* f_facets[4];
            int n_facets = 0;
            for (auto f_facet : frontFacets)
            {
                if (!*f_facet)
                    continue;

                f_facets[n_facets] = f_facet;
                n_facets++;
                if (n_facets == 4)
                {
                    innerSimps.push_back((new Simplex3(
                        **(*(*f_facet)->facet->edges[0])->vertexes[0],
                        **(*(*f_facet)->facet->edges[0])->vertexes[1],
                        **(*f_facet)->facet->findVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]),
                        **(*f_facets[0])->facet->findVertexNotIncludedInEdge(**Facet3::intersectAlongAnEdge(*(*f_facets[0])->facet, *(*f_facet)->facet))))->getPtrToUPtr());
                    break;
                }
            }

            for (auto &f_facet : frontFacets)
            {
                if (*f_facet)
                    delete f_facet->release();
            }
            for (auto &f_edge : frontEdges)
            {
                if (*f_edge)
                    delete f_edge->release();
            }
            return;
        }

        if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
            removePtrsToNullptr(frontEdges);
        if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
            removePtrsToNullptr(frontFacets);
    }

    polycr->outputData();
    if (!frontContainsOfOnly1Simplex3OrEmpty())
        throw std::logic_error("Error in function: Crystallite::processAngles");
}

void Crystallite3::processAngles_OLD(Polycrystal3* polycr)
{
    if (frontContainsOfOnly1Simplex3OrEmpty())
    {
        throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");

        for (auto &f_facet : frontFacets)
        {
            if (*f_facet)
                delete f_facet->release();
        }
        for (auto &f_edge : frontEdges)
        {
            if (*f_edge)
                delete f_edge->release();
        }
        return;
    }

    double max_excos = 2.0;
    for (unique_ptr<FrontEdge3>* curr_f_edge = currentFrontEdge(max_excos); ; curr_f_edge = currentFrontEdge(max_excos))
    {
        if (!curr_f_edge)
        {
            polycr->outputData();
            throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");
        }

        //double curr_f_edge_magn = (*curr_f_edge)->edge->magnitude();
        //if (curr_f_edge_magn < 0.2 * _preferredLength ||
        //    curr_f_edge_magn > 2.5 * _preferredLength)
        //{
        //    //for (int i = 0; i < 3; i++)
        //    //{
        //        smoothAroundFrontVertex((*curr_f_edge)->edge->vertexes[0]);
        //        smoothAroundFrontVertex((*curr_f_edge)->edge->vertexes[1]);
        //    //}
        //    for (auto f_edge : frontEdges)
        //    {
        //        if (!*f_edge)
        //            continue;
        //
        //        if ((*f_edge)->edge->contains(**(*curr_f_edge)->edge->vertexes[0]) ||
        //            (*f_edge)->edge->contains(**(*curr_f_edge)->edge->vertexes[1]))
        //        {
        //            (*f_edge)->needProcessing = true;
        //        }
        //    }
        //
        //    max_excos = 2.0;
        //    continue;
        //}
        //
        if ((*curr_f_edge)->getAngleExCos(this) > COS_DEG_100)
        {
            if (parallelFacetsCheck(curr_f_edge) ||
                edgeIntersectionCheck(curr_f_edge) ||
                facetsIntersectionCheck(curr_f_edge) ||
                anyVertexInsidePotentialSimplex3Check(curr_f_edge) ||
                frontSplitCheck(curr_f_edge))
            {
                //if ((*curr_f_edge)->getAngleExCos(frontFacets) < COS_DEG_100 &&
                //    !frontSplitCheck(curr_f_edge))
                //{
                //    Vec3 new_vert_pos;
                //    if (edgeIntersectionCheck(curr_f_edge) ||
                //        parallelFacetsCheck(curr_f_edge) ||
                //        !NewVertexPosition_OLD(curr_f_edge, new_vert_pos)) // Experemental
                //    {
                //        max_excos = (*curr_f_edge)->getAngleExCos(frontFacets);
                //        continue;
                //    }
                //
                //    ExhaustWithNewVertex_OLD(curr_f_edge, new_vert_pos);
                //    max_excos = 2.0;
                //}
                //else 
                //{
                //    max_excos = (*curr_f_edge)->getAngleExCos(frontFacets);
                //}
                //continue;
        
                max_excos = (*curr_f_edge)->getAngleExCos(this);
                continue;
            }
        
            exhaustWithoutNewVertex(curr_f_edge);
            max_excos = 2.0;
        }
        else
        {
            Vec3 new_vert_pos;
            if (edgeIntersectionCheck(curr_f_edge) ||
                parallelFacetsCheck(curr_f_edge) ||
                !NewVertexPosition_OLD(curr_f_edge, new_vert_pos))
            {
                max_excos = (*curr_f_edge)->getAngleExCos(this);
                continue;
            }
        
            ExhaustWithNewVertex_OLD(curr_f_edge, new_vert_pos);
            max_excos = 2.0;
        }
        
        //polycr->outputData();

        if (frontContainsOfOnly1Simplex3OrEmpty())
        {
            unique_ptr<FrontFacet3>* f_facets[4];
            int n_facets = 0;
            for (auto f_facet : frontFacets)
            {
                if (!*f_facet)
                    continue;

                f_facets[n_facets] = f_facet;
                n_facets++;
                if (n_facets == 4)
                {
                    innerSimps.push_back((new Simplex3(
                        **(*(*f_facet)->facet->edges[0])->vertexes[0],
                        **(*(*f_facet)->facet->edges[0])->vertexes[1],
                        **(*f_facet)->facet->findVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]),
                        **(*f_facets[0])->facet->findVertexNotIncludedInEdge(**Facet3::intersectAlongAnEdge(*(*f_facets[0])->facet, *(*f_facet)->facet))))->getPtrToUPtr());
                    break;
                }
            }

            for (auto &f_facet : frontFacets)
            {
                if (*f_facet)
                    delete f_facet->release();
            }
            for (auto &f_edge : frontEdges)
            {
                if (*f_edge)
                    delete f_edge->release();
            }
            return;
        }

        if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
            removePtrsToNullptr(frontEdges);
        if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
            removePtrsToNullptr(frontFacets);
    }

    polycr->outputData();
    if (!frontContainsOfOnly1Simplex3OrEmpty())
        throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");
}

void Crystallite3::generateMesh(const double preferredLength, Polycrystal3* polycr)
{
    _preferredLength = preferredLength;
    ////polycr->outputData();
    //if (ProcessVerySmallAngles(polycr))
    //{
    //    if (globalIntersectionCheck())
    //        std::cout << "Bad!";
    //    return;
    //}
    ////polycr->outputData();
    //if (ProcessSmallAngles(polycr))
    //{
    //    if (globalIntersectionCheck())
    //        std::cout << "Bad!";
    //    return;
    //}
    ////polycr->outputData();
    //if (ProcessMediumAngles(polycr))
    //{
    //    if (globalIntersectionCheck())
    //        std::cout << "Bad!";
    //    return;
    //}
    //polycr->outputData();
    //ProcessLargeAngles(polycr);
    computeFrontNormals();
    processAngles(polycr);
    //polycr->outputData();
    if (globalIntersectionCheck())
        throw std::logic_error("Intersection error. Crystallite3::globalIntersectionCheck returned true.");
    smoothMesh(10);
}

const bool Crystallite3::globalIntersectionCheck()
{
    for (auto &edge : innerEdges)
    {
        if (*edge &&
            edgeGlobalIntersectionCheck(edge))
            return true;
    }

    return false;
}

void Crystallite3::smoothMesh(int iterationsNum)
{
    for (int i = 0; i < iterationsNum; i++)
        for (auto &vert : innerVerts)
        {
            if (!*vert)
                continue;

            Vec3 shift;
            int delta_shifts_num = 0;
            for (auto &edge : innerEdges)
            {
                if (!*edge)
                    continue;

                if (vert == (*edge)->vertexes[0])
                {
                    shift += **(*edge)->vertexes[1] - **vert;
                }
                else if (vert == (*edge)->vertexes[1])
                {
                    shift += **(*edge)->vertexes[0] - **vert;
                }
                delta_shifts_num++;
            }
            shift /= delta_shifts_num;
            (*vert)->setPosition((*vert)->getPosition() + shift);
        }
}

void Crystallite3::smoothNotFinisedMesh(int iterationsNum)
{
    for (int i = 0; i < iterationsNum; i++)
        for (auto &vert : innerVerts)
        {
            if (!*vert)
                continue;

            smoothAroundFrontVertex(vert);
        }
}

void Crystallite3::smoothFront(int iterationsNum)
{
    for (int i = 0; i < iterationsNum; i++)
        for (auto &vert : innerVerts)
        {
            if (!*vert)
                continue;

            smoothAroundFrontVertex(vert);
        }
}

void Crystallite3::smoothAroundFrontVertex(unique_ptr<Vertex3>* frontVert)
{
    if ((*frontVert)->belongsToShellFacet ||
        (*frontVert)->belongsToShellEdge ||
        (*frontVert)->belongsToShellVertex)
    {
        return;
    }

    Vec3 shift;
    int delta_shifts_num = 0;
    for (auto &edge : innerEdges)
    {
        if (!*edge)
            continue;

        Vec3 d_shift;
        if (frontVert == (*edge)->vertexes[0])
        {
            d_shift = **(*edge)->vertexes[1] - **frontVert;
        }
        else if (frontVert == (*edge)->vertexes[1])
        {
            d_shift = **(*edge)->vertexes[0] - **frontVert;
        }
        shift += d_shift * (d_shift.magnitude() - _preferredLength);
        delta_shifts_num++;
    }
    shift /= delta_shifts_num;
    (*frontVert)->setPosition((*frontVert)->getPosition() + shift);
}

void Crystallite3::analyzeMeshQuality(double &out_minQuality, double &out_avQuality)
{
    size_t simps_num = 0ull;
    double av_q = 0.0;
    double min_q = 1.0;
    for (auto &simp : innerSimps)
    {
        if (!*simp)
            continue;

        double q = (*simp)->computeQuality();
        av_q += q;
        simps_num++;
        if (q < min_q)
            min_q = q;
    }
    av_q /= simps_num;

    out_minQuality = min_q;
    out_avQuality = av_q;
}

Crystallite3::Crystallite3() {}

Crystallite3::~Crystallite3()
{
    for (auto &simp : innerSimps)
    {
        if (*simp)
            delete simp->release();
        delete simp;
    }
    for (auto &facet : innerFacets)
    {
        if (*facet)
            delete facet->release();
        delete facet;
    }
    for (auto &edge : innerEdges)
    {
        if (*edge)
            delete edge->release();
        delete edge;
    }
    for (auto &vert : innerVerts)
    {
        if (*vert)
            delete vert->release();
        delete vert;
    }

    for (auto &ptr : frontFacets)
    {
        if (*ptr)
            delete ptr->release();
        delete ptr;
    }
    for (auto &ptr : frontEdges)
    {
        if (*ptr)
            delete ptr->release();
        delete ptr;
    }
}