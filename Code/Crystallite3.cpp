#include "Crystallite3.h"
#include <algorithm>
#include <iostream>
#include "SpatialAlgs.h"




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




bool Crystallite3::shellContainsVert(const Vertex3* vert)
{
    if (vert->belongsToShellFacet)
    {
        for (auto& s_facet : m_shellFacets)
        {
            if (s_facet == vert->belongsToShellFacet)
                return true;
        }
    }
    else if (vert->belongsToShellEdge)
    {
        for (auto& s_edge : m_shellEdges)
        {
            if (s_edge == vert->belongsToShellEdge)
                return true;
        }
    }
    else if (vert->belongsToShellVertex)
    {
        for (auto& s_edge : m_shellEdges)
        {
            if (s_edge->verts[0] == vert->belongsToShellVertex ||
                s_edge->verts[1] == vert->belongsToShellVertex)
                return true;
        }
    }

    return false;
}


void Crystallite3::setStartFront(const std::list<Edge3*>& edges, const std::list<Facet3*>& facets)
{
    for (auto& edge : edges)
    {
        if (shellContainsVert(edge->verts[0]) &&
            shellContainsVert(edge->verts[1]))
        {
            m_frontEdges.push_back(new FrontEdge3(this, edge));
        }
    }

    for (auto& facet : facets)
    {
        if (shellContainsVert(facet->edges[0]->verts[0]) &&
            shellContainsVert(facet->edges[0]->verts[1]) &&
            shellContainsVert(facet->findVertNotIncludedInEdge(facet->edges[0])))
        {
            m_frontFacets.push_back(new FrontFacet3(this, facet));
        }
    }
}


void Crystallite3::computeFrontNormals()
{
    for (auto& facet : m_frontFacets)
        facet->computeNormal();
}




ShellEdge3* Crystallite3::findShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const
{
    for (auto& s_edge : m_shellEdges)
    {
        if ((s_edge->verts[0] == v0  &&
             s_edge->verts[1] == v1) ||
            (s_edge->verts[1] == v0  &&
             s_edge->verts[0] == v1))
            return s_edge;
    }

    return nullptr;
}




FrontFacet3* Crystallite3::findFrontFacet(const Facet3* facet)
{
    for (auto& f_facet : m_frontFacets)
    {
        if (f_facet->facet == facet)
            return f_facet;
    }

    return nullptr;
}


FrontEdge3* Crystallite3::findFrontEdge(const Vertex3* v0, const Vertex3* v1)
{
    for (auto& f_edge : m_frontEdges)
    {
        if (((f_edge->edge->verts[0] == v0 &&
              f_edge->edge->verts[1] == v1) ||
             (f_edge->edge->verts[1] == v0 &&
              f_edge->edge->verts[0] == v1)))
            return f_edge;
    }

    return nullptr;
}


FrontEdge3* Crystallite3::findFrontEdge(const Edge3* edge)
{
    for (auto& f_edge : m_frontEdges)
    {
        if (f_edge->edge == edge)
            return f_edge;
    }

    return nullptr;
}




FrontFacet3* Crystallite3::addFrontFacet3(const Edge3* edge0, const Edge3* edge1, const Edge3* edge2)
{
    FrontFacet3* new_f_facet = new FrontFacet3(this, new Facet3(edge0, edge1, edge2));

    m_frontFacets.push_back(new_f_facet);
    m_innerFacets.push_back(new_f_facet->facet);

    return new_f_facet;
}


FrontEdge3* Crystallite3::addFrontEdge3(const Vertex3* vert0, const Vertex3* vert1)
{
    FrontEdge3* new_f_edge = new FrontEdge3(this, new Edge3(vert0, vert1));

    m_frontEdges.push_back(new_f_edge);
    m_innerEdges.push_back(new_f_edge->edge);

    return new_f_edge;
}


void Crystallite3::addFrontEdge3(const FrontEdge3* frontEdge)
{
    m_frontEdges.push_back((FrontEdge3*)frontEdge);
    m_innerEdges.push_back(frontEdge->edge);
}




bool Crystallite3::vertInsideFrontCheck(const Vec3& v)
{
    int inters_num = 0;
    for (auto& f_facet : m_frontFacets)
    {       
        Vec3 dir = 0.3333333333333333 * (f_facet->computeCenter() - 3.0 * v);

        for (auto& f_facetj : m_frontFacets)
        {
            if (tva::spatialalgs::isRayIntersectTriangle(
                    v, dir,
                    f_facetj->facet->edges[0]->verts[0]->getPos(),
                    f_facetj->facet->edges[0]->verts[1]->getPos(),
                    f_facetj->facet->findVertNotIncludedInEdge(f_facetj->facet->edges[0])->getPos()))
            {
                inters_num++;
            }
        }

        break;
    }

    return inters_num % 2 == 1;
}


bool Crystallite3::lineSegmentGlobalIntersectionCheck(const Vec3& v0, const Vec3& v1)
{
    for (auto& facet : m_innerFacets)
    {
        if (tva::spatialalgs::isSegmentIntersectTriangle(
                v0, v1,
                facet->edges[0]->verts[0]->getPos(),
                facet->edges[0]->verts[1]->getPos(),
                facet->findVertNotIncludedInEdge(facet->edges[0])->getPos()))
            return true;
    }

    return false;
}


bool Crystallite3::segmentFrontIntersectionCheck(const Vec3& v0, const Vec3& v1)
{
    for (auto& f_facet : m_frontFacets)
    {
        if (tva::spatialalgs::isSegmentIntersectTriangle(
                v0, v1,
                f_facet->facet->edges[0]->verts[0]->getPos(),
                f_facet->facet->edges[0]->verts[1]->getPos(),
                f_facet->facet->findVertNotIncludedInEdge(f_facet->facet->edges[0])->getPos()))
            return true;
    }

    return false;
}


bool Crystallite3::edgeGlobalIntersectionCheck(const Edge3* edge)
{
    Vec3 delta = 1e-3 * (*edge->verts[1] - *edge->verts[0]);
    Vec3 segment[2];
    segment[0] = edge->verts[0]->getPos() + delta;
    segment[1] = edge->verts[1]->getPos() - delta;

    for (auto& facet : m_innerFacets)
    {
        if (!facet->contains(edge) &&
            tva::spatialalgs::isSegmentIntersectTriangle(
                segment[0], segment[1],
                facet->edges[0]->verts[0]->getPos(), 
                facet->edges[0]->verts[1]->getPos(), 
                facet->findVertNotIncludedInEdge(facet->edges[0])->getPos()))
            return true;
    }

    return false;
}


bool XOR(bool b0, bool b1)
{
    return (b0 || b1) && !(b0 && b1);
}


bool Crystallite3::edgeIntersectionCheck(FrontEdge3* frontEdge)
{
    auto opp_verts = frontEdge->findOppVerts();
    Vec3 opp_verts_poses[2];
    opp_verts_poses[0] = opp_verts.first->getPos();
    opp_verts_poses[1] = opp_verts.second->getPos();
    Vec3 opp_vert0_to_1 = opp_verts_poses[1] - opp_verts_poses[0];
    Vec3 delta_vec_0th_to_1st = 1e-3 * opp_vert0_to_1;

    if (segmentFrontIntersectionCheck(
        opp_verts_poses[0] + delta_vec_0th_to_1st,
        opp_verts_poses[1] - delta_vec_0th_to_1st))
        return true;

    for (auto& f_edge : m_frontEdges)
    {
        if (f_edge->edge->contains(opp_verts.first) &&
            f_edge->edge->contains(opp_verts.second))
            return false;
    }

    for (auto& f_edge : m_frontEdges)
    {
        bool contains[2];
        contains[0] = f_edge->edge->contains(opp_verts.first);
        contains[1] = f_edge->edge->contains(opp_verts.second);

        Vertex3* vert_buf;
        if (contains[0])
        {
            if (f_edge->edge->verts[0] == opp_verts.first)
                vert_buf = f_edge->edge->verts[1];
            else
                vert_buf = f_edge->edge->verts[0];

            if (vert_buf->getPos().distanceToSegment(opp_verts_poses[0], opp_verts_poses[1]) < 4e-3 * m_preferredLength)
                return true;
        }
        else if (contains[1])
        {
            if (f_edge->edge->verts[0] == opp_verts.second)
                vert_buf = f_edge->edge->verts[1];
            else
                vert_buf = f_edge->edge->verts[0];

            if (vert_buf->getPos().distanceToSegment(opp_verts_poses[0], opp_verts_poses[1]) < 4e-3 * m_preferredLength)
                return true;
        }
        else
        {
            if (tva::spatialalgs::segmentsDistance(
                    opp_verts_poses[0], opp_verts_poses[1],
                    f_edge->edge->verts[0]->getPos(), f_edge->edge->verts[1]->getPos()) < 4e-3 * m_preferredLength)
                return true;
        }
    }

    return false;
}


bool Crystallite3::facetIntersectionCheck(const Vertex3* v0, const Vertex3* v1, const Vec3& v2)
{
    for (auto& f_edge : m_frontEdges)
    {
        bool contains[2];
        contains[0] = f_edge->edge->contains(v0);
        contains[1] = f_edge->edge->contains(v1);

        if (contains[0] || contains[1])
            continue;

        Vec3 first_to_second_delta = NOT_TOO_CLOSE * (f_edge->edge->verts[1]->getPos() - f_edge->edge->verts[0]->getPos());
        Vec3 f_edge_verts_poses[2];
        f_edge_verts_poses[0] = f_edge->edge->verts[0]->getPos() - first_to_second_delta;
        f_edge_verts_poses[1] = f_edge->edge->verts[1]->getPos() + first_to_second_delta;

        if (tva::spatialalgs::isSegmentIntersectTriangle(
            f_edge_verts_poses[0], f_edge_verts_poses[1],
            v0->getPos(), v1->getPos(), v2))
            return true;
    }

    return false;
}


bool Crystallite3::facetIntersectionCheck(const Vertex3* v0, const Vertex3* v1, const Vertex3* v2)
{
    Vec3 verts_poses[3];
    verts_poses[0] = v0->getPos();
    verts_poses[1] = v1->getPos();
    verts_poses[2] = v2->getPos();

    for (auto& f_edge : m_frontEdges)
    {
        bool contains[3];
        contains[0] = f_edge->edge->contains(v0);
        contains[1] = f_edge->edge->contains(v1);
        contains[2] = f_edge->edge->contains(v2);

        if ((contains[0] && contains[1]) ||
            (contains[0] && contains[2]) ||
            (contains[1] && contains[2]))
            continue;

        Vec3 first_to_second_delta = 1e-3 * (f_edge->edge->verts[1]->getPos() - f_edge->edge->verts[0]->getPos());
        Vec3 f_edge_verts_poses[2];
        if (contains[0])
        {
            if (f_edge->edge->verts[0] == v0)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->getPos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->getPos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->getPos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->getPos() - first_to_second_delta;
            }
        }
        else if (contains[1])
        {
            if (f_edge->edge->verts[0] == v1)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->getPos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->getPos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->getPos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->getPos() - first_to_second_delta;
            }
        }
        else if (contains[2])
        {
            if (f_edge->edge->verts[0] == v2)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->getPos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->getPos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->getPos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->getPos() - first_to_second_delta;
            }
        }
        else
        {
            f_edge_verts_poses[0] = f_edge->edge->verts[0]->getPos();
            f_edge_verts_poses[1] = f_edge->edge->verts[1]->getPos();
        }

        if (tva::spatialalgs::isSegmentIntersectTriangle(
                f_edge_verts_poses[0], f_edge_verts_poses[1],
                verts_poses[0], verts_poses[1], verts_poses[2]))
            return true;
    }

    return false;
}


bool Crystallite3::facetsIntersectionCheck(FrontEdge3* frontEdge)
{
    auto opp_verts = frontEdge->findOppVerts();
        
    if (facetIntersectionCheck(frontEdge->edge->verts[0], opp_verts.first, opp_verts.second) ||
        facetIntersectionCheck(frontEdge->edge->verts[1], opp_verts.first, opp_verts.second))
        return true;

    return false;
}


bool Crystallite3::insideSimplex3Check(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& vert)
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

    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4] * 1.000001;
}


bool Crystallite3::anyVertInsidePotentialSimp3Check(FrontEdge3* frontEdge)
{
    auto opp_verts = frontEdge->findOppVerts();

    Vec3 points[4];
    points[0] = opp_verts.first->getPos();
    points[1] = opp_verts.second->getPos();
    points[2] = frontEdge->edge->verts[0]->getPos();
    points[3] = frontEdge->edge->verts[1]->getPos();
    for (auto& vert : m_innerVerts)
    {
        if (vert != opp_verts.first &&
            vert != opp_verts.second &&
            vert != frontEdge->edge->verts[0] &&
            vert != frontEdge->edge->verts[1] &&
            insideSimplex3Check(points[0], points[1], points[2], points[3], vert->getPos()))
            return true;
    }

    return false;
}


bool Crystallite3::frontSplitCheck(FrontEdge3* frontEdge)
{
    FrontEdge3* opp_f_edge = frontEdge->findOppEdge();
    if (!opp_f_edge)
        return false;

    auto opp_f_edge_opp_verts = opp_f_edge->findOppVerts();

    if (opp_f_edge_opp_verts.first == frontEdge->edge->verts[0] ||
        opp_f_edge_opp_verts.first == frontEdge->edge->verts[1] ||
        opp_f_edge_opp_verts.second == frontEdge->edge->verts[0] ||
        opp_f_edge_opp_verts.second == frontEdge->edge->verts[1])
        return false;

    return true;
}


bool Crystallite3::parallelFacetsCheck(FrontEdge3* frontEdge) const
{
    auto adj_f_facets = frontEdge->getAdjFacets();

    Vertex3* opp_verts[2];
    opp_verts[0] = adj_f_facets.first->facet->findVertNotIncludedInEdge(frontEdge->edge);
    opp_verts[1] = adj_f_facets.second->facet->findVertNotIncludedInEdge(frontEdge->edge);

    Vec3 plane0[2];
    plane0[0] = *opp_verts[0] - *frontEdge->edge->verts[0];
    plane0[1] = *opp_verts[1] - *frontEdge->edge->verts[0];
    Vec3 normal0 = Vec3::crossProduct(plane0[0], plane0[1]).normalize();
    Vec3 plane1[2];
    plane1[0] = *opp_verts[0] - *frontEdge->edge->verts[1];
    plane1[1] = *opp_verts[1] - *frontEdge->edge->verts[1];
    Vec3 normal1 = Vec3::crossProduct(plane1[0], plane1[1]).normalize();

    for (auto& f_facet : m_frontFacets)
    {
        Edge3* inters_reses[2];
        if ((f_facet != adj_f_facets.first) &&
            (f_facet != adj_f_facets.second))
        {
            inters_reses[0] = Facet3::intersectAlongAnEdge(f_facet->facet, adj_f_facets.first->facet);
            inters_reses[1] = Facet3::intersectAlongAnEdge(f_facet->facet, adj_f_facets.second->facet);
            if (!XOR((bool)inters_reses[0], (bool)inters_reses[1]))
                continue;

            Vertex3* f_facet_to_verts[2];
            f_facet_to_verts[0] = f_facet->facet->edges[0]->verts[0];
            f_facet_to_verts[1] = f_facet->facet->edges[0]->verts[1];
            Vertex3* f_facet_from_vert = f_facet->facet->findVertNotIncludedInEdge(f_facet->facet->edges[0]);

            Vec3 f_plane[2];
            f_plane[0] = *f_facet_to_verts[0] - *f_facet_from_vert;
            f_plane[1] = *f_facet_to_verts[1] - *f_facet_from_vert;
            Vec3 f_normal = Vec3::crossProduct(f_plane[0], f_plane[1]).normalize();

            if (abs(abs(Vec3::dotProduct(f_normal, normal0)) - 1.0) < 1e-6 ||
                abs(abs(Vec3::dotProduct(f_normal, normal1)) - 1.0) < 1e-6)
            {
                int i = inters_reses[0] ? 0 : 1;

                Vec3 border_verts[2];
                border_verts[0] = inters_reses[i]->verts[0]->getPos();
                border_verts[1] = inters_reses[i]->verts[1]->getPos();

                Vec3 main_facet_3rd_vert;
                if (inters_reses[i]->contains(opp_verts[0]))
                    main_facet_3rd_vert = opp_verts[1]->getPos();
                else
                    main_facet_3rd_vert = opp_verts[0]->getPos();

                Vec3 curr_facet_3rd_vert = f_facet->facet->findVertNotIncludedInEdge(inters_reses[i])->getPos();

                Vec3 main_facet_cross = Vec3::crossProduct(main_facet_3rd_vert - border_verts[0], main_facet_3rd_vert - border_verts[1]);
                Vec3 curr_facet_cross = Vec3::crossProduct(curr_facet_3rd_vert - border_verts[0], curr_facet_3rd_vert - border_verts[1]);
                if (Vec3::dotProduct(main_facet_cross, curr_facet_cross) > 0.0)
                    return true;
            }
        }
    }

    return false;
}




double Crystallite3::computeMinEdgesLength(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
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


double Crystallite3::computeMaxEdgesLength(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
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


double Crystallite3::computeSimp3SimpleQuality(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return sqrt(sqr_min_max.first / sqr_min_max.second);
}


double Crystallite3::computeSimp3SimpleSqrQuality(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return sqr_min_max.first / sqr_min_max.second;
}




FrontEdge3* Crystallite3::currentFrontEdge(double maxExCos) const
{
    double curr_max_excos = -2.0;
    FrontEdge3* curr_max_f_edge = nullptr;
    for (auto& f_edge : m_frontEdges)
    {
        double curr_excos = f_edge->getAngleExCos();
        if (curr_excos > curr_max_excos &&
            curr_excos < maxExCos)
        {
            curr_max_excos = curr_excos;
            curr_max_f_edge = f_edge;
        }
    }

    return curr_max_f_edge;
}


bool Crystallite3::exhaustWithoutNewVertPriorityPredicate(FrontEdge3* currentFrontEdge)
{
    if (currentFrontEdge->getAngleExCos() > COS_DEG_70)
        return true;

    auto opp_verts = currentFrontEdge->findOppVerts();
    if (findFrontEdge(opp_verts.first, opp_verts.second) ||
        (currentFrontEdge->getAngleExCos() < COS_DEG_70 &&
         currentFrontEdge->getAngleExCos() > COS_DEG_100 &&
         (*opp_verts.second - *opp_verts.first).sqrMagnitude() <= m_preferredLength * m_preferredLength))
        return true;

    return false;
}


bool Crystallite3::exhaustWithNewVertPriorityPredicate(FrontEdge3* currentFrontEdge)
{
    if (currentFrontEdge->getAngleExCos() < COS_DEG_120)
        return true;

    return false;
}

 
Crystallite3::ExhaustType Crystallite3::computeExhaustionTypeQualityPriority(
    FrontEdge3* currentFrontEdge, 
    FrontFacet3** out_withNWFrontFacet, Vec3** out_withNWNewVertPos)
{
    if (frontSplitCheck(currentFrontEdge))
        return DONT_EXHAUST;
    
    if (parallelFacetsCheck(currentFrontEdge) ||
        edgeIntersectionCheck(currentFrontEdge) ||
        facetsIntersectionCheck(currentFrontEdge) ||
        anyVertInsidePotentialSimp3Check(currentFrontEdge))
    {
        return WITH_NEW_VERTEX;
    }

    auto opp_verts = currentFrontEdge->findOppVerts();
    double without_nv_quality = computeSimp3SimpleSqrQuality(
        currentFrontEdge->edge->verts[0]->getPos(),
        currentFrontEdge->edge->verts[1]->getPos(),
        opp_verts.first->getPos(),
        opp_verts.second->getPos());

    FrontFacet3* f_facet = chooseFacetForExhaustionWithNewVert(currentFrontEdge);
    Vec3 new_vert_pos;
    if (!tryComputeNewVertPos(f_facet, new_vert_pos))
        return DONT_EXHAUST;

    double with_nv_quality = computeSimp3SimpleSqrQuality(
        f_facet->facet->edges[0]->verts[0]->getPos(),
        f_facet->facet->edges[0]->verts[1]->getPos(),
        f_facet->facet->findVertNotIncludedInEdge(f_facet->facet->edges[0])->getPos(),
        new_vert_pos);

    if (without_nv_quality > with_nv_quality)
        return WITHOUT_NEW_VERTEX;

    *out_withNWFrontFacet = f_facet;
    *out_withNWNewVertPos = new Vec3(new_vert_pos);
    return WITH_NEW_VERTEX;
}




Vec3 Crystallite3::computeNormalInSimp3(FrontFacet3* frontFacet, const Vec3& oppositeVertPos)
{
    Vec3 third_pos = frontFacet->facet->findVertNotIncludedInEdge(frontFacet->facet->edges[0])->getPos();
    Vec3 normal = Vec3::crossProduct(
        frontFacet->facet->edges[0]->verts[0]->getPos() - third_pos,
        frontFacet->facet->edges[0]->verts[1]->getPos() - third_pos).normalize();
    if (Vec3::dotProduct(normal, oppositeVertPos - third_pos) > 0.0)
        normal *= -1.0;

    return normal;
}


Vec3 Crystallite3::computeNormalInSimp3(FrontFacet3* frontFacet, Edge3* oneOfRemainingEdges)
{
    Vec3 opposite = frontFacet->facet->contains(oneOfRemainingEdges->verts[0]) ?
        oneOfRemainingEdges->verts[1]->getPos() :
        oneOfRemainingEdges->verts[0]->getPos();
    Vec3 third_pos = frontFacet->facet->findVertNotIncludedInEdge(frontFacet->facet->edges[0])->getPos();
    Vec3 normal = Vec3::crossProduct(
        frontFacet->facet->edges[0]->verts[0]->getPos() - third_pos,
        frontFacet->facet->edges[0]->verts[1]->getPos() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - third_pos) > 0.0)
        normal *= -1.0;

    return normal;
}




void Crystallite3::exhaustWithoutNewVertOppEdgeExists(FrontEdge3* frontEdge, FrontEdge3* oppositeEdge)
{
    FrontEdge3* opp_f_edge = oppositeEdge;

    auto opp_f_facets = opp_f_edge->getAdjFacets();

    FrontFacet3* main_f_facets[3];
    Vertex3* main_vert;
    auto f_edge_adj_facets = frontEdge->getAdjFacets();
    main_f_facets[0] = f_edge_adj_facets.first;
    main_f_facets[1] = f_edge_adj_facets.second;

    if (opp_f_facets.first->facet->contains(frontEdge->edge->verts[0]))
    {
        main_f_facets[2] = opp_f_facets.first;
        main_vert = frontEdge->edge->verts[0];
    }
    else if (opp_f_facets.first->facet->contains(frontEdge->edge->verts[1]))
    {
        main_f_facets[2] = opp_f_facets.first;
        main_vert = frontEdge->edge->verts[1];
    }
    else if (opp_f_facets.second->facet->contains(frontEdge->edge->verts[0]))
    {
        main_f_facets[2] = opp_f_facets.second;
        main_vert = frontEdge->edge->verts[0];
    }
    else
    {
        main_f_facets[2] = opp_f_facets.second;
        main_vert = frontEdge->edge->verts[1];
    }

    FrontEdge3* new_f_facet_f_edge;
    Edge3* new_simp_edges[3];
    for (int i = 0; i < 3; i++)
    {
        new_simp_edges[i] = main_f_facets[i]->facet->findEdgeNotContainingVert(main_vert);
        new_f_facet_f_edge = findFrontEdge(new_simp_edges[i]);
        new_f_facet_f_edge->refreshAngleData();
        new_f_facet_f_edge->removeAdjFacetFromPair(main_f_facets[i]);
    }

    auto new_f_facet = addFrontFacet3(
        new_simp_edges[0],
        new_simp_edges[1],
        new_simp_edges[2]);
    Vec3 center = new_f_facet->computeCenter();
    Vec3 opposite = main_vert->getPos();
    Vec3 third_pos = new_f_facet->facet->findVertNotIncludedInEdge(oppositeEdge->edge)->getPos();
    Vec3 normal = Vec3::crossProduct(
        oppositeEdge->edge->verts[0]->getPos() - third_pos,
        oppositeEdge->edge->verts[1]->getPos() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    new_f_facet->setNormal(normal);

    m_innerSimps.push_back(new Simplex3(
        frontEdge->edge->verts[0],
        frontEdge->edge->verts[1],
        opp_f_edge->edge->verts[0],
        opp_f_edge->edge->verts[1]));

    FrontEdge3* f_edge_ = nullptr;
    for (auto& f_facet : main_f_facets)
    {
        for (auto& edge : f_facet->facet->edges)
        {
            if (edge->contains(main_vert) &&
                (f_edge_ = findFrontEdge(edge)))
            {
                m_frontEdges.erase(std::find(m_frontEdges.begin(), m_frontEdges.end(), f_edge_));
                delete f_edge_;
            }
        }
    }

    for (auto& f_facet : main_f_facets)
    {
        m_frontFacets.erase(std::find(m_frontFacets.begin(), m_frontFacets.end(), f_facet));
        delete f_facet;
    }
}

//#pragma optimize( "", off )
void Crystallite3::exhaustWithoutNewVertOppEdgeDontExists(FrontEdge3* frontEdge)
{
    // This volatile helps to avoid computational error.
    volatile auto adj_f_facets = frontEdge->getAdjFacets();

    Vertex3* opp_verts[2];
    opp_verts[0] = adj_f_facets.first->facet->findVertNotIncludedInEdge(frontEdge->edge);
    opp_verts[1] = adj_f_facets.second->facet->findVertNotIncludedInEdge(frontEdge->edge);
    
    FrontEdge3* opp_f_edge = addFrontEdge3(opp_verts[0], opp_verts[1]);

    Edge3* new_simp_edges[3];
    new_simp_edges[2] = opp_f_edge->edge;

    FrontEdge3* buf_f_edge;
    for (auto& edge : adj_f_facets.first->facet->edges)
    {
        if (edge->contains(frontEdge->edge->verts[0]) &&
            (edge->contains(opp_f_edge->edge->verts[0]) ||
             edge->contains(opp_f_edge->edge->verts[1])))
        {
            new_simp_edges[0] = edge;
            buf_f_edge = findFrontEdge(edge);
            buf_f_edge->refreshAngleData();
            buf_f_edge->removeAdjFacetFromPair(adj_f_facets.first);
            buf_f_edge->removeAdjFacetFromPair(adj_f_facets.second);
            break;
        }
    }
    for (auto& edge : adj_f_facets.second->facet->edges)
    {
        if (edge->contains(frontEdge->edge->verts[0]) &&
            (edge->contains(opp_f_edge->edge->verts[0]) ||
             edge->contains(opp_f_edge->edge->verts[1])))
        {
            new_simp_edges[1] = edge;
            buf_f_edge = findFrontEdge(edge);
            buf_f_edge->refreshAngleData();
            buf_f_edge->removeAdjFacetFromPair(adj_f_facets.first);
            buf_f_edge->removeAdjFacetFromPair(adj_f_facets.second);
            break;
        }
    }
    auto new_f_facet = addFrontFacet3(
        new_simp_edges[0],
        new_simp_edges[1],
        new_simp_edges[2]);
    Vec3 center = new_f_facet->computeCenter();
    Vec3 opposite = new_simp_edges[0]->contains(frontEdge->edge->verts[0]) ?
        frontEdge->edge->verts[1]->getPos() :
        frontEdge->edge->verts[0]->getPos();
    Vec3 third_pos = new_f_facet->facet->findVertNotIncludedInEdge(opp_f_edge->edge)->getPos();
    Vec3 normal = Vec3::crossProduct(
        opp_f_edge->edge->verts[0]->getPos() - third_pos,
        opp_f_edge->edge->verts[1]->getPos() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    new_f_facet->setNormal(normal);

    for (auto& edge : adj_f_facets.first->facet->edges)
    {
        if (edge->contains(frontEdge->edge->verts[1]) &&
            (edge->contains(opp_f_edge->edge->verts[0]) ||
             edge->contains(opp_f_edge->edge->verts[1])))
        {
            new_simp_edges[0] = edge;
            buf_f_edge = findFrontEdge(edge);
            buf_f_edge->refreshAngleData();
            buf_f_edge->removeAdjFacetFromPair(adj_f_facets.first);
            buf_f_edge->removeAdjFacetFromPair(adj_f_facets.second);
            break;
        }
    }
    for (auto& edge : adj_f_facets.second->facet->edges)
    {
        if (edge->contains(frontEdge->edge->verts[1]) &&
            (edge->contains(opp_f_edge->edge->verts[0]) ||
             edge->contains(opp_f_edge->edge->verts[1])))
        {
            new_simp_edges[1] = edge;
            buf_f_edge = findFrontEdge(edge);
            buf_f_edge->refreshAngleData();
            buf_f_edge->removeAdjFacetFromPair(adj_f_facets.first);
            buf_f_edge->removeAdjFacetFromPair(adj_f_facets.second);
            break;
        }
    }
    new_f_facet = addFrontFacet3(
        new_simp_edges[0],
        new_simp_edges[1],
        new_simp_edges[2]);
    center = new_f_facet->computeCenter();
    opposite = new_simp_edges[0]->contains(frontEdge->edge->verts[0]) ?
        frontEdge->edge->verts[1]->getPos() :
        frontEdge->edge->verts[0]->getPos();
    third_pos = new_f_facet->facet->findVertNotIncludedInEdge(opp_f_edge->edge)->getPos();
    normal = Vec3::crossProduct(
        opp_f_edge->edge->verts[0]->getPos() - third_pos,
        opp_f_edge->edge->verts[1]->getPos() - third_pos).normalize();
    if (Vec3::dotProduct(normal, opposite - center) > 0.0)
        normal *= -1.0;
    new_f_facet->setNormal(normal);

    m_innerSimps.push_back(new Simplex3(
        frontEdge->edge->verts[0],
        frontEdge->edge->verts[1],
        opp_verts[0],
        opp_verts[1]));

    
    m_frontFacets.erase(std::find(m_frontFacets.begin(), m_frontFacets.end(), adj_f_facets.first));
    m_frontFacets.erase(std::find(m_frontFacets.begin(), m_frontFacets.end(), adj_f_facets.second));
    delete adj_f_facets.first;
    delete adj_f_facets.second;

    m_frontEdges.erase(std::find(m_frontEdges.begin(), m_frontEdges.end(), frontEdge));
    delete frontEdge;
}
//#pragma optimize( "", on )

void Crystallite3::exhaustWithoutNewVert(FrontEdge3* frontEdge, bool oppositeEdgeExistence, FrontEdge3* oppositeEdge)
{
    FrontEdge3* opp_f_edge = nullptr;
    if (oppositeEdgeExistence && oppositeEdge)
        opp_f_edge = oppositeEdge;
    else if (oppositeEdgeExistence)
        opp_f_edge = frontEdge->findOppEdge();

    if ((oppositeEdge && oppositeEdgeExistence) ||
        opp_f_edge)
    {
        exhaustWithoutNewVertOppEdgeExists(frontEdge, opp_f_edge);
        //std::cout << "__00__";
    }
    else
    {
        exhaustWithoutNewVertOppEdgeDontExists(frontEdge);
        //std::cout << "__01__";
    }
}




bool Crystallite3::tryComputeNewVertPosType3(FrontFacet3* frontFacet, Vec3& out_pos)
{
    Edge3* main_edges[2];
    main_edges[0] = frontFacet->facet->edges[0];
    main_edges[1] = frontFacet->facet->edges[1];
    FrontEdge3* main_f_edges[2];
    main_f_edges[0] = findFrontEdge(main_edges[0]);
    main_f_edges[1] = findFrontEdge(main_edges[1]);
    auto main_vert = main_edges[0]->verts[0] == main_edges[1]->verts[0] ?
        main_edges[0]->verts[0] :
        (main_edges[0]->verts[0] == main_edges[1]->verts[1] ?
            main_edges[0]->verts[0] :
            main_edges[0]->verts[1]);
    auto third_edge = frontFacet->facet->findEdgeNotContainingVert(main_vert);
    auto third_f_edge = findFrontEdge(third_edge);

    double av_magn = 0.3333333333333333 * (
          frontFacet->facet->edges[0]->magnitude()
        + frontFacet->facet->edges[1]->magnitude()
        + frontFacet->facet->edges[2]->magnitude());

    auto v0 = main_edges[0]->verts[0] == main_vert ? main_edges[0]->verts[1] : main_edges[0]->verts[0];
    auto v1 = main_edges[1]->verts[0] == main_vert ? main_edges[1]->verts[1] : main_edges[1]->verts[0];
    auto v2 = main_vert;

    auto adj_f_facets = main_f_edges[0]->getAdjFacets();
    auto fn0 = adj_f_facets.first == frontFacet ? adj_f_facets.second : adj_f_facets.first;
    adj_f_facets = main_f_edges[1]->getAdjFacets();
    auto fn1 = adj_f_facets.first == frontFacet ? adj_f_facets.second : adj_f_facets.first;
    adj_f_facets = third_f_edge->getAdjFacets();
    auto fn2 = adj_f_facets.first == frontFacet ? adj_f_facets.second : adj_f_facets.first;

    const Vec3& v0_pos = v0->getPos();
    const Vec3& v1_pos = v1->getPos();
    const Vec3& v2_pos = v2->getPos();

    Vec3 n_m = frontFacet->getNormal();
    Vec3 n_n0 = fn0->getNormal();
    Vec3 n_n1 = fn1->getNormal();
    Vec3 n_n2 = fn2->getNormal();

    Vec3 e_mn0 = (n_m + n_n0).normalize();
    Vec3 e_mn1 = (n_m + n_n1).normalize();
    Vec3 e_mn2 = n_m + n_n2;

    Vec3 np_mn0 = Vec3::crossProduct(v2_pos - v0_pos, e_mn0);
    Vec3 np_mn1 = Vec3::crossProduct(v2_pos - v1_pos, e_mn1);

    Vec3 e = Vec3::crossProduct(np_mn0, np_mn1);

    Vec3 new_pos = tva::spatialalgs::lineIntersectPlane(v2_pos, e, v0_pos, v1_pos, v0_pos + e_mn2);

    Vec3 v0_to_np = new_pos - v0_pos;
    Vec3 v1_to_np = new_pos - v1_pos;
    Vec3 v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + 1e-3 * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + 1e-3 * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + 1e-3 * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
        !vertInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

    out_pos = new_pos;
    return true;
}


bool Crystallite3::tryComputeNewVertPosType2(const int smallAngleIndex0, const double angleCos0, const int smallAngleIndex1, const double angleCos1, FrontFacet3* frontFacet, Vec3& out_pos)
{
    Edge3* main_edges[2];
    main_edges[0] = frontFacet->facet->edges[smallAngleIndex0];
    main_edges[1] = frontFacet->facet->edges[smallAngleIndex1];
    FrontEdge3* main_f_edges[2];
    main_f_edges[0] = findFrontEdge(main_edges[0]);
    main_f_edges[1] = findFrontEdge(main_edges[1]);
    auto main_vert = main_edges[0]->verts[0] == main_edges[1]->verts[0] ?
        main_edges[0]->verts[0] :
        (main_edges[0]->verts[0] == main_edges[1]->verts[1] ?
            main_edges[0]->verts[0] :
            main_edges[0]->verts[1]);

    double av_magn = 0.3333333333333333 * (
          frontFacet->facet->edges[0]->magnitude()
        + frontFacet->facet->edges[1]->magnitude()
        + frontFacet->facet->edges[2]->magnitude());

    auto v0 = main_edges[0]->verts[0] == main_vert ? main_edges[0]->verts[1] : main_edges[0]->verts[0];
    auto v1 = main_edges[1]->verts[0] == main_vert ? main_edges[1]->verts[1] : main_edges[1]->verts[0];
    auto v2 = main_vert;

    auto adj_f_facets = main_f_edges[0]->getAdjFacets();
    auto fn0 = adj_f_facets.first == frontFacet ? adj_f_facets.second : adj_f_facets.first;
    adj_f_facets = main_f_edges[1]->getAdjFacets();
    auto fn1 = adj_f_facets.first == frontFacet ? adj_f_facets.second : adj_f_facets.first;

    const Vec3& v0_pos = v0->getPos();
    const Vec3& v1_pos = v1->getPos();
    const Vec3& v2_pos = v2->getPos();

    Vec3 n_m = frontFacet->getNormal();
    Vec3 n_n0 = fn0->getNormal();
    Vec3 n_n1 = fn1->getNormal();
    
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
        !vertInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

    out_pos = new_pos;
    return true;
}


bool Crystallite3::tryComputeNewVertPosType1(const int smallAngleIndex, const double angleCos, FrontFacet3* frontFacet, Vec3& out_pos)
{
    auto main_edge = frontFacet->facet->edges[smallAngleIndex];
    auto main_f_edge = findFrontEdge(main_edge);
    Vec3 main_edge_poses[2];
    main_edge_poses[0] = main_edge->verts[0]->getPos();
    main_edge_poses[1] = main_edge->verts[1]->getPos();

    double av_magn = 0.3333333333333333 * (
          frontFacet->facet->edges[0]->magnitude()
        + frontFacet->facet->edges[1]->magnitude()
        + frontFacet->facet->edges[2]->magnitude());

    auto v0 = main_edge->verts[0];
    auto v1 = main_edge->verts[1];
    auto v2 = frontFacet->facet->findVertNotIncludedInEdge(main_edge);

    auto adj_f_facets = main_f_edge->getAdjFacets();
    auto fn = adj_f_facets.first == frontFacet ? adj_f_facets.second : adj_f_facets.first;
    auto vn = fn->facet->findVertNotIncludedInEdge(main_edge);

    const Vec3& v0_pos = v0->getPos();
    const Vec3& v1_pos = v1->getPos();
    const Vec3& v2_pos = v2->getPos();
    const Vec3& vn_pos = vn->getPos();

    Vec3 v2pr = tva::spatialalgs::project(v2_pos, v0_pos, v1_pos);
    Vec3 vnpr = tva::spatialalgs::project(vn_pos, v0_pos, v1_pos);

    Vec3 c = 0.25 * (v0_pos + v1_pos + v2pr + vnpr);
    Vec3 e = (frontFacet->getNormal() + fn->getNormal()).normalize();
    Vec3 new_pos = angleCos > COS_DEG_80 ?
        c + 0.6 * ONE_PLUS_SQRT2_DIV_SQRT3 * av_magn * e :
        c + 0.9 * ONE_PLUS_SQRT2_DIV_SQRT3 * av_magn * e;

    Vec3 v0_to_np = new_pos - v0_pos;
    Vec3 v1_to_np = new_pos - v1_pos;
    Vec3 v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + 1e-3 * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + 1e-3 * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + 1e-3 * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
        !vertInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

    out_pos = new_pos;
    return true;
}


bool Crystallite3::tryComputeNewVertPosType0(FrontFacet3* frontFacet, Vec3& out_pos)
{
    double av_magn = 0.3333333333333333 * (
          frontFacet->facet->edges[0]->magnitude()
        + frontFacet->facet->edges[1]->magnitude()
        + frontFacet->facet->edges[2]->magnitude());
    Vec3 new_pos = frontFacet->computeCenter() + SQRT_2_DIV_3 * av_magn * frontFacet->getNormal();

    auto v0 = frontFacet->facet->edges[0]->verts[0];
    auto v1 = frontFacet->facet->edges[0]->verts[1];
    auto v2 = frontFacet->facet->findVertNotIncludedInEdge(frontFacet->facet->edges[0]);
    const Vec3& v0_pos = v0->getPos();
    const Vec3& v1_pos = v1->getPos();
    const Vec3& v2_pos = v2->getPos();
    Vec3 v0_to_np = new_pos - v0_pos;
    Vec3 v1_to_np = new_pos - v1_pos;
    Vec3 v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + 1e-3 * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + 1e-3 * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + 1e-3 * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
        !vertInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

    out_pos = new_pos;
    return true;
}


bool Crystallite3::tryComputeNewVertPos(FrontFacet3* frontFacet, Vec3& out_pos)
{
    double angs_coses[3]
    { 
        findFrontEdge(frontFacet->facet->edges[0])->getAngleExCos(),
        findFrontEdge(frontFacet->facet->edges[1])->getAngleExCos(),
        findFrontEdge(frontFacet->facet->edges[2])->getAngleExCos()
    };
    int indexes[3];
    int small_angs_num = 0;
    if (angs_coses[0] > COS_DEG_140) indexes[small_angs_num++] = 0;
    if (angs_coses[1] > COS_DEG_140) indexes[small_angs_num++] = 1;
    if (angs_coses[2] > COS_DEG_140) indexes[small_angs_num++] = 2;

    switch (small_angs_num)
    {
    case 0: return tryComputeNewVertPosType0(frontFacet, out_pos);        
    case 1: return tryComputeNewVertPosType1(indexes[0], angs_coses[indexes[0]], frontFacet, out_pos);
    case 2: return tryComputeNewVertPosType2(indexes[0], angs_coses[indexes[0]], indexes[1], angs_coses[indexes[1]], frontFacet, out_pos);
    case 3: return tryComputeNewVertPosType3(frontFacet, out_pos);
    }

    return true;
}




FrontFacet3* Crystallite3::chooseFacetForExhaustionWithNewVert(FrontEdge3* frontEdge)
{
    auto adj_f_facets = frontEdge->getAdjFacets();

    return adj_f_facets.first->computeQuality() < adj_f_facets.second->computeQuality() ?
        adj_f_facets.first :
        adj_f_facets.second;
}


void Crystallite3::exhaustWithNewVert(FrontFacet3* frontFacet, Vec3 vertPos)
{
    Vertex3* new_vert = new Vertex3(vertPos);
    m_innerVerts.push_back(new_vert);

    Edge3* new_simp_edges[6];
    new_simp_edges[0] = frontFacet->facet->edges[0];
    auto far_vert = frontFacet->facet->findVertNotIncludedInEdge(new_simp_edges[0]);
    new_simp_edges[1] = frontFacet->facet->findEdge(new_simp_edges[0]->verts[1], far_vert);
    new_simp_edges[2] = frontFacet->facet->findEdge(new_simp_edges[0]->verts[0], far_vert);
    new_simp_edges[3] = addFrontEdge3(new_simp_edges[0]->verts[0], new_vert)->edge;
    new_simp_edges[4] = addFrontEdge3(new_simp_edges[0]->verts[1], new_vert)->edge;
    new_simp_edges[5] = addFrontEdge3(far_vert, new_vert)->edge;

    FrontEdge3* new_simp_f_edges[3];
    new_simp_f_edges[0] = findFrontEdge(new_simp_edges[0]);
    new_simp_f_edges[1] = findFrontEdge(new_simp_edges[1]);
    new_simp_f_edges[2] = findFrontEdge(new_simp_edges[2]);

    new_simp_f_edges[0]->refreshAngleData();
    new_simp_f_edges[1]->refreshAngleData();
    new_simp_f_edges[2]->refreshAngleData();

    new_simp_f_edges[0]->removeAdjFacetFromPair(frontFacet);
    new_simp_f_edges[1]->removeAdjFacetFromPair(frontFacet);
    new_simp_f_edges[2]->removeAdjFacetFromPair(frontFacet);
    
    auto new_f_facet = addFrontFacet3(
        new_simp_edges[0],
        new_simp_edges[3],
        new_simp_edges[4]);
    new_f_facet->setNormal(computeNormalInSimp3(new_f_facet, far_vert->getPos()));

    new_f_facet = addFrontFacet3(
        new_simp_edges[2],
        new_simp_edges[3],
        new_simp_edges[5]);
    new_f_facet->setNormal(computeNormalInSimp3(new_f_facet, new_simp_edges[0]->verts[1]->getPos()));

    new_f_facet = addFrontFacet3(
        new_simp_edges[1],
        new_simp_edges[4],
        new_simp_edges[5]);
    new_f_facet->setNormal(computeNormalInSimp3(new_f_facet, new_simp_edges[0]->verts[0]->getPos()));

    m_innerSimps.push_back(new Simplex3(
        new_simp_edges[0]->verts[0],
        new_simp_edges[0]->verts[1],
        new_simp_edges[5]->verts[0],
        new_simp_edges[5]->verts[1]));

    m_frontFacets.erase(std::find(m_frontFacets.begin(), m_frontFacets.end(), frontFacet));
    delete frontFacet;
}




bool Crystallite3::tryExhaustWithoutNewVert(FrontEdge3* frontEdge, bool oppositeEdgeExistence, FrontEdge3* oppositeEdge)
{
    if (parallelFacetsCheck(frontEdge) ||
        edgeIntersectionCheck(frontEdge) ||
        facetsIntersectionCheck(frontEdge) ||
        anyVertInsidePotentialSimp3Check(frontEdge) ||
        frontSplitCheck(frontEdge))
    {
        return false;
    }

    exhaustWithoutNewVert(frontEdge, oppositeEdgeExistence, oppositeEdge);
    return true;
}


bool Crystallite3::tryExhaustWithNewVert(FrontEdge3* frontEdge)
{
    if (edgeIntersectionCheck(frontEdge) ||
        parallelFacetsCheck(frontEdge))
    {
        return false;
    }

    auto exhaust_f_facet = chooseFacetForExhaustionWithNewVert(frontEdge);
    Vec3 new_vert_pos;
    if (!tryComputeNewVertPos(exhaust_f_facet, new_vert_pos))
        return false;

    exhaustWithNewVert(exhaust_f_facet, new_vert_pos);
    return true;
}




void Crystallite3::processAngles()
{
    if (m_frontFacets.size() <= 4ull)
    {
        throw std::logic_error("Error in function: Crystallite::processAngles");

        for (auto& f_facet : m_frontFacets)
            delete f_facet;
        m_frontFacets.clear();

        for (auto& f_edge : m_frontEdges)
            delete f_edge;
        m_frontEdges.clear();

        return;
    }

    double max_excos = 2.0;
    for (FrontEdge3* curr_f_edge = currentFrontEdge(max_excos);; curr_f_edge = currentFrontEdge(max_excos))
    {
        if (!curr_f_edge)
            throw std::logic_error("Crystallite3::currentFrontEdge returned nullptr");
        
        if (exhaustWithoutNewVertPriorityPredicate(curr_f_edge))
        {
            if (!tryExhaustWithoutNewVert(curr_f_edge))
            {
                max_excos = curr_f_edge->getAngleExCos();
                continue;
            }
        }
        else if (exhaustWithNewVertPriorityPredicate(curr_f_edge))
        {
            if (!tryExhaustWithNewVert(curr_f_edge))
            {
                max_excos = curr_f_edge->getAngleExCos();
                continue;
            }
        }
        else
        {
            FrontFacet3* exhaust_from_f_facet = nullptr;
            Vec3* new_vert_pos = nullptr;
            switch (computeExhaustionTypeQualityPriority(curr_f_edge, &exhaust_from_f_facet, &new_vert_pos))
            {
            case WITHOUT_NEW_VERTEX:
                exhaustWithoutNewVert(curr_f_edge);
                break;

            case WITH_NEW_VERTEX:
                if (new_vert_pos)
                {
                    exhaustWithNewVert(exhaust_from_f_facet, *new_vert_pos);
                    delete new_vert_pos;
                }
                else
                {
                    if (!tryExhaustWithNewVert(curr_f_edge))
                    {
                        max_excos = curr_f_edge->getAngleExCos();
                        continue;
                    }
                }
                break;

            case DONT_EXHAUST:
                max_excos = curr_f_edge->getAngleExCos();
                continue;
                break;
            }
        }
        max_excos = 2.0;

        if (m_frontFacets.size() == 4ull)
        {
            FrontFacet3* f_facets[4];
            int n_facets = 0;
            for (auto& f_facet : m_frontFacets)
            {
                f_facets[n_facets] = f_facet;
                n_facets++;
                if (n_facets == 4)
                {
                    m_innerSimps.push_back(new Simplex3(
                        f_facet->facet->edges[0]->verts[0],
                        f_facet->facet->edges[0]->verts[1],
                        f_facet->facet->findVertNotIncludedInEdge(f_facet->facet->edges[0]),
                        f_facets[0]->facet->findVertNotIncludedInEdge(Facet3::intersectAlongAnEdge(f_facets[0]->facet, f_facet->facet))));
                    break;
                }
            }

            for (auto& f_facet : m_frontFacets)
                delete f_facet;
            m_frontFacets.clear();

            for (auto& f_edge : m_frontEdges)
                delete f_edge;
            m_frontEdges.clear();

            return;
        }
    }
}


void Crystallite3::generateMesh(const double preferredLength)
{
    m_preferredLength = preferredLength;
    computeFrontNormals();
    processAngles();
    if (globalIntersectionCheck())
        throw std::logic_error("Intersection error. Crystallite3::globalIntersectionCheck returned true.");
    smoothMesh(10);
}


bool Crystallite3::globalIntersectionCheck()
{
    for (auto& edge : m_innerEdges)
    {
        if (edgeGlobalIntersectionCheck(edge))
            return true;
    }

    return false;
}


void Crystallite3::smoothMesh(int iterationsNum)
{
    for (int i = 0; i < iterationsNum; i++)
        for (auto &vert : m_innerVerts)
        {
            Vec3 shift;
            int delta_shifts_num = 0;
            for (auto &edge : m_innerEdges)
            {
                if (vert == edge->verts[0])
                {
                    shift += *edge->verts[1] - *vert;
                }
                else if (vert == edge->verts[1])
                {
                    shift += *edge->verts[0] - *vert;
                }
                delta_shifts_num++;
            }
            shift /= delta_shifts_num;
            vert->setPos(vert->getPos() + shift);
        }
}


void Crystallite3::smoothNotFinisedMesh(int iterationsNum)
{
    for (int i = 0; i < iterationsNum; i++)
        for (auto &vert : m_innerVerts)
            smoothAroundFrontVert(vert);
}


void Crystallite3::smoothFront(int iterationsNum)
{
    for (int i = 0; i < iterationsNum; i++)
        for (auto &vert : m_innerVerts)
            smoothAroundFrontVert(vert);
}


void Crystallite3::smoothAroundFrontVert(Vertex3* fVert)
{
    if (fVert->belongsToShellFacet ||
        fVert->belongsToShellEdge ||
        fVert->belongsToShellVertex)
    {
        return;
    }

    Vec3 shift;
    int delta_shifts_num = 0;
    for (auto &edge : m_innerEdges)
    {
        Vec3 d_shift;
        if (fVert == edge->verts[0])
        {
            d_shift = *edge->verts[1] - *fVert;
        }
        else if (fVert == edge->verts[1])
        {
            d_shift = *edge->verts[0] - *fVert;
        }
        shift += d_shift * (d_shift.magnitude() - m_preferredLength);
        delta_shifts_num++;
    }
    shift /= delta_shifts_num;
    fVert->setPos(fVert->getPos() + shift);
}


std::pair<double, double> Crystallite3::analyzeMeshQuality()
{
    size_t simps_num = 0ull;
    double av_q = 0.0;
    double min_q = 1.0;
    for (auto &simp : m_innerSimps)
    {
        double q = simp->computeQuality();
        av_q += q;
        simps_num++;
        if (q < min_q)
            min_q = q;
    }
    av_q /= simps_num;

    return { min_q, av_q };
}




Crystallite3::Crystallite3() {}


Crystallite3::~Crystallite3()
{
    for (auto& simp : m_innerSimps)
        delete simp;
    for (auto& facet : m_innerFacets)
        delete facet;
    for (auto& edge : m_innerEdges)
        delete edge;
    for (auto& vert : m_innerVerts)
        delete vert;
}



void Crystallite3::addShellFacet(const ShellFacet3* shellFacet)
{
    m_shellFacets.push_back((ShellFacet3*)shellFacet);
}

void Crystallite3::addShellEdge(const ShellEdge3* shellEdge)
{
    m_shellEdges.push_back((ShellEdge3*)shellEdge);
}

bool Crystallite3::shellFacetsContains(const ShellFacet3* shellFacet)
{
    return std::find(m_shellFacets.begin(), m_shellFacets.end(), shellFacet) != m_shellFacets.end();
}


bool Crystallite3::shellEdgesContains(const ShellEdge3* shellEdge)
{
    return std::find(m_shellEdges.begin(), m_shellEdges.end(), shellEdge) != m_shellEdges.end();
}

const std::vector<Simplex3*>& Crystallite3::getInnerSimplexes3()
{
    return m_innerSimps;
}

const std::vector<Facet3*>& Crystallite3::getInnerFacets()
{
    return m_innerFacets;
}

const std::vector<Vertex3*>& Crystallite3::getInnerVertexes()
{
    return m_innerVerts;
}

const std::list<FrontFacet3*>& Crystallite3::getFrontFacets()
{
    return m_frontFacets;
}

const std::list<FrontEdge3*>& Crystallite3::getFrontEdges()
{
    return m_frontEdges;
}
