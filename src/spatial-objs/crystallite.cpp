#include "crystallite.h"
#include <algorithm>
#include <iostream>
#include <float.h>
#include "helpers/spatial-algs/spatial-algs.h"

#include "helpers/cosd-values.h"

using namespace pmg;
using tva::Vec;
using tva::Point;

using FrSuFacet = front::surface::Facet;
using FrSuEdge  = front::surface::Edge;
using pair_ff   = std::pair<FrSuFacet*, FrSuFacet*>;


#define DET(a, b, c, d) \
        (a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
        (std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
        (BETWEEN(corner0[0], corner1[0], point[0]) && \
         BETWEEN(corner0[1], corner1[1], point[1]))

#define ALPHA_P      70.52877936550931
#define DEG_1_IN_RAD  0.0174532925199432957

#define SQRT3_2              0.8660254037844386
#define SQRT_2_3             0.8164965809277260
#define ONE_PLUS_SQRT2_SQRT3 1.3938468501173517

#define NOT_TOO_CLOSE          3e-1
#define FROM_VERT_COEF         1e-2
#define EDGES_INTERS_DIST_COEF 4e-3

#define K_MAXD 0.2
#define K_D    0.3


template <typename T>
constexpr double degToRad(T value)
{
    return value * DEG_1_IN_RAD;
}


//#define DEV_DEBUG




bool Crystallite::shellContains(const Vertex* vert) const
{
    if (vert->belongsToShellFacet)
    {
        for (auto& s_facet : m_shellFacets)
            if (s_facet == vert->belongsToShellFacet)
                return true;
    }
    else if (vert->belongsToShellEdge)
    {
        for (auto& s_edge : m_shellEdges)
            if (s_edge == vert->belongsToShellEdge)
                return true;
    }
    else if (vert->belongsToShellVertex)
    {
        for (auto& s_edge : m_shellEdges)
            if (s_edge->verts[0] == vert->belongsToShellVertex ||
                s_edge->verts[1] == vert->belongsToShellVertex)
                return true;
    }

    return false;
}




void Crystallite::initializeFFacetFEdges(FrSuFacet* fFacet) const
{
    int n_added = 0;
    for (auto& fedge : m_frontEdges)
    {
        if (fFacet->facet->contains(fedge->edge))
        {
            fFacet->addFEdge(fedge);
            if (++n_added == 3)
                return;
        }
    }

    throw std::logic_error("Crystallite::initializeFFacetFEdges didn't find 3 front edges.");
}


void Crystallite::initializeFront()
{
    for (auto& sedge : m_shellEdges)
        for (auto& edge : sedge->innerEdges())
            m_frontEdges.push_back(new FrSuEdge(this, edge));

    for (auto& sfacet : m_shellFacets)
        for (auto& edge : sfacet->innerEdges())
            m_frontEdges.push_back(new FrSuEdge(this, edge));

    for (auto& sfacet : m_shellFacets)
        for (auto& facet : sfacet->innerFacets())
            m_frontFacets.push_back(new FrSuFacet(this, facet));

    for (auto& ffacet : m_frontFacets)
        initializeFFacetFEdges(ffacet);
}


void Crystallite::computeFrontNormals()
{
    for (auto& facet : m_frontFacets)
        facet->computeNormal();
}




shell::Edge* Crystallite::findShellEdge(const shell::Vertex* v0, const shell::Vertex* v1) const
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


FrSuFacet* Crystallite::findFrontFacet(const Facet* facet) const
{
    for (auto& f_facet : m_frontFacets)
    {
        if (f_facet->facet == facet)
            return f_facet;
    }

    return nullptr;
}


std::vector<FrSuEdge*> Crystallite::findFEdge(const Vertex* v0, const Vertex* v1) const
{
    std::vector<FrSuEdge*> res;
    for (auto& f_edge : m_frontEdges)
    {
        if (((f_edge->edge->verts[0] == v0 &&
              f_edge->edge->verts[1] == v1) ||
             (f_edge->edge->verts[1] == v0 &&
              f_edge->edge->verts[0] == v1)))
        {
            res.push_back(f_edge);
        }
    }

    return res;
}


std::vector<FrSuEdge*> Crystallite::findFEdge(const Edge* edge) const
{
    std::vector<FrSuEdge*> res;
    for (auto& f_edge : m_frontEdges)
    {
        if (f_edge->edge == edge)
        {
            res.push_back(f_edge);
        }
    }

    return res;
}




FrSuFacet* Crystallite::addToFront(const Facet* facet, bool addInner)
{
    FrSuFacet* new_ffacet = new FrSuFacet(this, facet);
    m_frontFacets.push_back(new_ffacet);
    if (addInner)
        m_innerFacets.push_back(new_ffacet->facet);
    return new_ffacet;
}


FrSuEdge* Crystallite::addToFront(const pmg::Edge* edge, bool addInner)
{
    FrSuEdge* new_f_edge = new FrSuEdge(this, edge);
    m_frontEdges.push_back(new_f_edge);
    if (addInner)
        m_innerEdges.push_back(new_f_edge->edge);
    return new_f_edge;
}




void Crystallite::removeFromFront(FrSuFacet* fFacet)
{
    m_frontFacets.erase(std::find(m_frontFacets.begin(), m_frontFacets.end(), fFacet));
    delete fFacet;
}


void Crystallite::removeFromFront(FrSuEdge* fEdge)
{
    m_frontEdges.erase(std::find(m_frontEdges.begin(), m_frontEdges.end(), fEdge));
    delete fEdge;
}




bool Crystallite::vertInsideFrontCheck(const Vec& v) const
{
    int inters_num = 0;
    Vec dir = 0.3333333333333333 * (m_frontFacets.front()->computeCenter() - 3.0 * v);

    for (auto& f_facet : m_frontFacets)
    {
        if (tva::spatalgs::doesRayIntersectTriangle(
                v, dir,
                f_facet->facet->edges[0]->verts[0]->pos(),
                f_facet->facet->edges[0]->verts[1]->pos(),
                f_facet->facet->findVertNot(f_facet->facet->edges[0])->pos()))
        {
             inters_num++;
        }
    }

    return inters_num % 2 == 1;
}


bool Crystallite::segmentGlobalIntersectionCheck(const Vec& v0, const Vec& v1) const
{
    for (auto& facet : m_innerFacets)
    {
        if (tva::spatalgs::doesSegmentIntersectTriangle(
                v0, v1,
                facet->edges[0]->verts[0]->pos(),
                facet->edges[0]->verts[1]->pos(),
                facet->findVertNot(facet->edges[0])->pos()))
            return true;
    }

    return false;
}


bool Crystallite::segmentFrontIntersectionCheck(const Vec& v0, const Vec& v1) const
{
    for (auto& f_facet : m_frontFacets)
    {
        if (tva::spatalgs::doesSegmentIntersectTriangle(
                v0, v1,
                f_facet->facet->edges[0]->verts[0]->pos(),
                f_facet->facet->edges[0]->verts[1]->pos(),
                f_facet->facet->findVertNot(f_facet->facet->edges[0])->pos()))
            return true;
    }

    return false;
}


bool Crystallite::edgeGlobalIntersectionCheck(const Edge* edge) const
{
    Vec delta = 8.0 * FROM_VERT_COEF * (*edge->verts[1] - *edge->verts[0]);
    Vec segment[2];
    segment[0] = edge->verts[0]->pos() + delta;
    segment[1] = edge->verts[1]->pos() - delta;

    for (auto& facet : m_innerFacets)
    {
        if (!facet->contains(edge) &&
            tva::spatalgs::doesSegmentIntersectTriangle(
                segment[0], segment[1],
                facet->edges[0]->verts[0]->pos(),
                facet->edges[0]->verts[1]->pos(),
                facet->findVertNot(facet->edges[0])->pos()))
            return true;
    }

    return false;
}


bool XOR(bool b0, bool b1)
{
    return (b0 || b1) && !(b0 && b1);
}


bool Crystallite::edgeIntersectionCheck(FrSuEdge* frontEdge) const
{
    auto opp_verts = frontEdge->findOppVerts();
    Vec opp_verts_poses[2];
    opp_verts_poses[0] = std::get<0>(opp_verts)->pos();
    opp_verts_poses[1] = std::get<1>(opp_verts)->pos();
    Vec opp_vert0_to_1 = opp_verts_poses[1] - opp_verts_poses[0];
    Vec delta_vec_0th_to_1st = FROM_VERT_COEF * opp_vert0_to_1;

    if (segmentFrontIntersectionCheck(
        opp_verts_poses[0] + delta_vec_0th_to_1st,
        opp_verts_poses[1] - delta_vec_0th_to_1st))
        return true;

    for (auto& f_edge : m_frontEdges)
    {
        if (f_edge->edge->contains(std::get<0>(opp_verts)) &&
            f_edge->edge->contains(std::get<1>(opp_verts)))
            return false;
    }

    for (auto& f_edge : m_frontEdges)
    {
        Vertex* vert_buf;
        if (bool contains[2] { f_edge->edge->contains(std::get<0>(opp_verts)),
                               f_edge->edge->contains(std::get<1>(opp_verts)) };
            contains[0])
        {
            if (f_edge->edge->verts[0] == std::get<0>(opp_verts))
                vert_buf = f_edge->edge->verts[1];
            else
                vert_buf = f_edge->edge->verts[0];

            if (tva::spatalgs::distancePointToSegment(vert_buf->pos(), opp_verts_poses[0], opp_verts_poses[1]) < EDGES_INTERS_DIST_COEF * m_preferredLength)
                return true;
        }
        else if (contains[1])
        {
            if (f_edge->edge->verts[0] == std::get<1>(opp_verts))
                vert_buf = f_edge->edge->verts[1];
            else
                vert_buf = f_edge->edge->verts[0];

            if (tva::spatalgs::distancePointToSegment(vert_buf->pos(), opp_verts_poses[0], opp_verts_poses[1]) < EDGES_INTERS_DIST_COEF * m_preferredLength)
                return true;
        }
        else
        {
            if (tva::spatalgs::segmentsDistance(
                    opp_verts_poses[0], opp_verts_poses[1],
                    f_edge->edge->verts[0]->pos(), f_edge->edge->verts[1]->pos()) < EDGES_INTERS_DIST_COEF * m_preferredLength)
                return true;
        }
    }

    return false;
}


bool Crystallite::facetIntersectionCheck(const Vertex* v0, const Vertex* v1, const Vec& v2) const
{
    for (auto& f_edge : m_frontEdges)
    {
        bool contains[2];
        contains[0] = f_edge->edge->contains(v0);
        contains[1] = f_edge->edge->contains(v1);

        if (contains[0] || contains[1])
            continue;

        Vec first_to_second_delta = NOT_TOO_CLOSE * (f_edge->edge->verts[1]->pos() - f_edge->edge->verts[0]->pos());
        Vec f_edge_verts_poses[2];
        f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos() - first_to_second_delta;
        f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos() + first_to_second_delta;

        if (tva::spatalgs::doesSegmentIntersectTriangle(
            f_edge_verts_poses[0], f_edge_verts_poses[1],
            v0->pos(), v1->pos(), v2))
            return true;
    }

    return false;
}


bool Crystallite::facetIntersectionCheck(const Vertex* v0, const Vertex* v1, const Vertex* v2) const
{
    Vec verts_poses[3];
    verts_poses[0] = v0->pos();
    verts_poses[1] = v1->pos();
    verts_poses[2] = v2->pos();

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

        Vec first_to_second_delta = FROM_VERT_COEF * (f_edge->edge->verts[1]->pos() - f_edge->edge->verts[0]->pos());
        Vec f_edge_verts_poses[2];
        if (contains[0])
        {
            if (f_edge->edge->verts[0] == v0)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos() - first_to_second_delta;
            }
        }
        else if (contains[1])
        {
            if (f_edge->edge->verts[0] == v1)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos() - first_to_second_delta;
            }
        }
        else if (contains[2])
        {
            if (f_edge->edge->verts[0] == v2)
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos() + first_to_second_delta;
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos();
            }
            else
            {
                f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos();
                f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos() - first_to_second_delta;
            }
        }
        else
        {
            f_edge_verts_poses[0] = f_edge->edge->verts[0]->pos();
            f_edge_verts_poses[1] = f_edge->edge->verts[1]->pos();
        }

        if (tva::spatalgs::doesSegmentIntersectTriangle(
                f_edge_verts_poses[0], f_edge_verts_poses[1],
                verts_poses[0], verts_poses[1], verts_poses[2]))
            return true;
    }

    return false;
}


bool Crystallite::facetsIntersectionCheck(FrSuEdge* frontEdge) const
{
    auto opp_verts = frontEdge->findOppVerts();
        
    if (facetIntersectionCheck(frontEdge->edge->verts[0], std::get<0>(opp_verts), std::get<1>(opp_verts)) ||
        facetIntersectionCheck(frontEdge->edge->verts[1], std::get<0>(opp_verts), std::get<1>(opp_verts)))
        return true;

    return false;
}


bool Crystallite::insideTetrCheck(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3, const Vec& vert) const
{
    Vec vert_to_p0 = p0 - vert;
    Vec vert_to_p1 = p1 - vert;
    Vec vert_to_p2 = p2 - vert;
    Vec vert_to_p3 = p3 - vert;

    double abs_mixed_prods[5];
    abs_mixed_prods[0] = std::abs(Vec::mixed(vert_to_p0, vert_to_p2, vert_to_p3));
    abs_mixed_prods[1] = std::abs(Vec::mixed(vert_to_p0, vert_to_p1, vert_to_p2));
    abs_mixed_prods[2] = std::abs(Vec::mixed(vert_to_p0, vert_to_p1, vert_to_p3));
    abs_mixed_prods[3] = std::abs(Vec::mixed(vert_to_p1, vert_to_p2, vert_to_p3));
    abs_mixed_prods[4] = std::abs(Vec::mixed(p1 - p0, p2 - p0, p3 - p0));

    return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4]/* * 1.000001*/;
}


bool Crystallite::anyVertInsidePotentialTetrCheck(FrSuEdge* fEdge) const
{
    auto opp_verts = fEdge->findOppVerts();

    Vec points[4];
    points[0] = std::get<0>(opp_verts)->pos();
    points[1] = std::get<1>(opp_verts)->pos();
    points[2] = fEdge->edge->verts[0]->pos();
    points[3] = fEdge->edge->verts[1]->pos();
    for (auto& vert : m_innerVerts)
    {
        if (vert != std::get<0>(opp_verts) &&
            vert != std::get<1>(opp_verts) &&
            vert != fEdge->edge->verts[0] &&
            vert != fEdge->edge->verts[1] &&
            insideTetrCheck(points[0], points[1], points[2], points[3], vert->pos()))
            return true;
    }

    return false;
}


bool Crystallite::frontSplitCheck(FrSuEdge* fEdge, FrSuEdge* oppFEdge) const
{
    auto opp_fedge = oppFEdge;
    if (!opp_fedge)
        opp_fedge = fEdge->findOppEdge();
    if (!opp_fedge)
        return false;

    auto opp_fedge_opp_verts = opp_fedge->findOppVerts();

    if (std::get<0>(opp_fedge_opp_verts) == fEdge->edge->verts[0] ||
        std::get<0>(opp_fedge_opp_verts) == fEdge->edge->verts[1] ||
        std::get<1>(opp_fedge_opp_verts) == fEdge->edge->verts[0] ||
        std::get<1>(opp_fedge_opp_verts) == fEdge->edge->verts[1])
        return false;

    return true;
}


bool Crystallite::frontCollapseCheck(FrSuEdge* fEdge, FrSuEdge* oppFEdge) const
{
    auto opp_fedge = oppFEdge;
    if (!opp_fedge)
        opp_fedge = fEdge->findOppEdge();
    if (!opp_fedge)
        return false;

    auto opp_ffacets = opp_fedge->getAdjFFacets();

    if (!fEdge->edge->contains(opp_ffacets.first->facet->findVertNot(opp_fedge->edge)) ||
        !fEdge->edge->contains(opp_ffacets.second->facet->findVertNot(opp_fedge->edge)))
        return false;

    return true;
}


bool Crystallite::parallelFacetsCheck(FrSuEdge* fEdge) const
{
    auto adj_f_facets = fEdge->getAdjFFacets();

    Vertex* opp_verts[2];
    opp_verts[0] = std::get<0>(adj_f_facets)->facet->findVertNot(fEdge->edge);
    opp_verts[1] = std::get<1>(adj_f_facets)->facet->findVertNot(fEdge->edge);

    Vec plane0[2];
    plane0[0] = *opp_verts[0] - *fEdge->edge->verts[0];
    plane0[1] = *opp_verts[1] - *fEdge->edge->verts[0];
    Vec normal0 = Vec::cross(plane0[0], plane0[1]).normalize();
    Vec plane1[2];
    plane1[0] = *opp_verts[0] - *fEdge->edge->verts[1];
    plane1[1] = *opp_verts[1] - *fEdge->edge->verts[1];
    Vec normal1 = Vec::cross(plane1[0], plane1[1]).normalize();

    for (auto& f_facet : m_frontFacets)
    {
        Edge* inters_reses[2];
        if ((f_facet != std::get<0>(adj_f_facets)) &&
            (f_facet != std::get<1>(adj_f_facets)))
        {
            inters_reses[0] = Facet::intersectAlongEdge(f_facet->facet, std::get<0>(adj_f_facets)->facet);
            inters_reses[1] = Facet::intersectAlongEdge(f_facet->facet, std::get<1>(adj_f_facets)->facet);
            if (!XOR(static_cast<bool>(inters_reses[0]), static_cast<bool>(inters_reses[1])))
                continue;

            Vertex* f_facet_to_verts[2];
            f_facet_to_verts[0] = f_facet->facet->edges[0]->verts[0];
            f_facet_to_verts[1] = f_facet->facet->edges[0]->verts[1];
            Vertex* f_facet_from_vert = f_facet->facet->findVertNot(f_facet->facet->edges[0]);

            Vec f_plane[2];
            f_plane[0] = *f_facet_to_verts[0] - *f_facet_from_vert;
            f_plane[1] = *f_facet_to_verts[1] - *f_facet_from_vert;
            Vec f_normal = Vec::cross(f_plane[0], f_plane[1]).normalize();

            if (std::abs(std::abs(Vec::dot(f_normal, normal0)) - 1.0) < 1e-6 ||
                std::abs(std::abs(Vec::dot(f_normal, normal1)) - 1.0) < 1e-6)
            {
                int i = inters_reses[0] ? 0 : 1;

                Vec border_verts[2];
                border_verts[0] = inters_reses[i]->verts[0]->pos();
                border_verts[1] = inters_reses[i]->verts[1]->pos();

                Vec main_facet_3rd_vert;
                if (inters_reses[i]->contains(opp_verts[0]))
                    main_facet_3rd_vert = opp_verts[1]->pos();
                else
                    main_facet_3rd_vert = opp_verts[0]->pos();

                Vec curr_facet_3rd_vert = f_facet->facet->findVertNot(inters_reses[i])->pos();

                Vec main_facet_cross = Vec::cross(main_facet_3rd_vert - border_verts[0], main_facet_3rd_vert - border_verts[1]);
                Vec curr_facet_cross = Vec::cross(curr_facet_3rd_vert - border_verts[0], curr_facet_3rd_vert - border_verts[1]);
                if (Vec::dot(main_facet_cross, curr_facet_cross) > 0.0)
                    return true;
            }
        }
    }

    return false;
}


bool Crystallite::doesFrontIntersectSphere(const Point& center, double radius) const
{
    for (auto& ffacet : m_frontFacets)
    {
        Point triangle[3];
        triangle[0] = ffacet->facet->edges[0]->verts[0]->pos();
        triangle[1] = ffacet->facet->edges[0]->verts[1]->pos();
        triangle[2] = ffacet->facet->findVertNot(ffacet->facet->edges[0])->pos();
        if (tva::spatalgs::doesTriangleIntersectSphere(triangle[0], triangle[1], triangle[2], center, radius))
            return true;
    }

    return false;
}




std::pair<double, double> Crystallite::computeMinMaxEdgesLengths(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3)
{
    auto min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    min_max.first = sqrt(min_max.first);
    min_max.second = sqrt(min_max.second);
    return min_max;
}


std::pair<double, double> Crystallite::computeMinMaxEdgesSqrLengths(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3)
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


double Crystallite::computeTetrSimpleQuality(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return sqrt(sqr_min_max.first / sqr_min_max.second);
}


double Crystallite::computeTetrSimpleSqrQuality(const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3)
{
    auto sqr_min_max = computeMinMaxEdgesSqrLengths(p0, p1, p2, p3);
    return sqr_min_max.first / sqr_min_max.second;
}




FrSuEdge* Crystallite::currentFrontEdge(double maxCompl) const
{
    double cur_max_compl = 0.0;
    FrSuEdge* cur_max_f_edge = nullptr;
    for (auto& f_edge : m_frontEdges)
    {
        double cur_compl = f_edge->complexity();
        if (cur_compl > cur_max_compl &&
            cur_compl < maxCompl)
        {
            cur_max_compl = cur_compl;
            cur_max_f_edge = f_edge;
        }
    }

    return cur_max_f_edge;
}


bool Crystallite::exhaustWithoutNewVertPriorityPredicate(FrSuEdge* curFEdge)
{
    if (curFEdge->angleExCos() > cosDeg<70>)
        return true;

    auto opp_verts = curFEdge->findOppVerts();
    if (curFEdge->findOppEdge() ||
        (curFEdge->angleExCos() < cosDeg<70> &&
         curFEdge->angleExCos() > cosDeg<100> &&
         (*std::get<1>(opp_verts) - *std::get<0>(opp_verts)).sqrMagnitude() <= m_preferredLength * m_preferredLength))
        return true;

    if (frontSplitCheck(curFEdge))
        return true;

    return false;
}


bool Crystallite::exhaustWithNewVertPriorityPredicate(FrSuEdge* currentFrontEdge)
{
    if (currentFrontEdge->angleExCos() < cosDeg<120>)
        return true;

    return false;
}

 
Crystallite::ExhaustType Crystallite::computeExhaustionTypeQualityPriority(
    FrSuEdge* currentFrontEdge,
    FrSuFacet*& out_withNWFrontFacet, Vec*& out_withNWNewVertPos)
{
    if (frontSplitCheck(currentFrontEdge))
        return ExhaustType::WithoutNewVert;
    
    if (parallelFacetsCheck(currentFrontEdge) ||
        edgeIntersectionCheck(currentFrontEdge) ||
        facetsIntersectionCheck(currentFrontEdge) ||
        anyVertInsidePotentialTetrCheck(currentFrontEdge))
    {
        return ExhaustType::WithNewVert;
    }

    auto opp_verts = currentFrontEdge->findOppVerts();
    double without_nv_quality = computeTetrSimpleSqrQuality(
        currentFrontEdge->edge->verts[0]->pos(),
        currentFrontEdge->edge->verts[1]->pos(),
        std::get<0>(opp_verts)->pos(),
        std::get<1>(opp_verts)->pos());

    FrSuFacet* f_facet = chooseFacetForExhaustionWithNewVert(currentFrontEdge);
    Vec new_vert_pos;
    if (!tryComputeNewVertPos(f_facet, new_vert_pos))
        return ExhaustType::DontExhaust;

    double with_nv_quality = computeTetrSimpleSqrQuality(
        f_facet->facet->edges[0]->verts[0]->pos(),
        f_facet->facet->edges[0]->verts[1]->pos(),
        f_facet->facet->findVertNot(f_facet->facet->edges[0])->pos(),
        new_vert_pos);

    if (without_nv_quality > with_nv_quality)
        return ExhaustType::WithoutNewVert;

    out_withNWFrontFacet = f_facet;
    out_withNWNewVertPos = new Vec(new_vert_pos);
    return ExhaustType::WithNewVert;
}




Vec Crystallite::computeNormalInTetr(const FrSuFacet* fFacet, const Vec& oppVertPos) const
{
    return computeNormalInTetr(fFacet->facet->edges[0]->verts[0]->pos(),
                               fFacet->facet->edges[0]->verts[1]->pos(),
                               fFacet->facet->findVertNot(fFacet->facet->edges[0])->pos(),
                               oppVertPos);
}


Vec Crystallite::computeNormalInTetr(const FrSuFacet* fFacet, const pmg::Edge* oneOfRemainingEdges) const
{
    Vec opp_pos = fFacet->facet->contains(oneOfRemainingEdges->verts[0]) ?
        oneOfRemainingEdges->verts[1]->pos() :
        oneOfRemainingEdges->verts[0]->pos();

    return computeNormalInTetr(fFacet, opp_pos);
}


Vec Crystallite::computeNormalInTetr(const Point& fFacetPos0, const Point& fFacetPos1, const Point& fFacetPos2, const Point& oppVertPos) const
{
    Vec normal = Vec::cross(
        fFacetPos0 - fFacetPos2,
        fFacetPos1 - fFacetPos2).normalize();
    if (Vec::dot(normal, oppVertPos - fFacetPos2) > 0.0)
        normal *= -1.0;

    return normal;
}




void Crystallite::setFEdgesInFrontSplit(const FrSuEdge* fEdge, FrSuEdge* newOppFEdges[2], FrSuFacet* newFFacets[2], pair_ff oppFFacets) const
{
//    double cpa_times[2][2];
//    cpa_times[0][0] = tva::spatalgs::cpaTime(newFFacets[0]->computeCenter(), newFFacets[0]->normal,
//                                             oppFFacets.first->computeCenter(), oppFFacets.first->normal);
//    cpa_times[0][1] = tva::spatalgs::cpaTime(newFFacets[0]->computeCenter(), newFFacets[0]->normal,
//                                             oppFFacets.second->computeCenter(), oppFFacets.second->normal);
//    cpa_times[1][0] = tva::spatalgs::cpaTime(newFFacets[1]->computeCenter(), newFFacets[1]->normal,
//                                             oppFFacets.first->computeCenter(), oppFFacets.first->normal);
//    cpa_times[1][1] = tva::spatalgs::cpaTime(newFFacets[1]->computeCenter(), newFFacets[1]->normal,
//                                             oppFFacets.second->computeCenter(), oppFFacets.second->normal);

//    std::pair<FrSuEdge*, FrSuFacet*> pair0;
//    std::pair<FrSuEdge*, FrSuFacet*> pair1;
//    bool trivial_solution = false;

//    if (cpa_times[0][0] > 1e-6 && cpa_times[1][1] > 1e-6)
//    {
//        pair0 = { newOppFEdges[0], oppFFacets.first };
//        pair1 = { newOppFEdges[1], oppFFacets.second };
//        trivial_solution = true;
//    }

//    if (cpa_times[0][1] > 1e-6 && cpa_times[1][0] > 1e-6)
//    {
//        pair0 = { newOppFEdges[0], oppFFacets.second };
//        pair1 = { newOppFEdges[1], oppFFacets.first };
//        trivial_solution = true;
//    }

//    if (trivial_solution)
//    {
//        pair0.first->fillAdjFFacets(newFFacets[0], pair0.second);
//        pair0.second->addFEdge(pair0.first);

//        pair1.first->fillAdjFFacets(newFFacets[1], pair1.second);
//        pair1.second->addFEdge(pair1.first);
//        return;
//    }

//    std::cout << "\nCrystallite::setFEdgesInFrontSplit here is not trivial solution...";
//    std::cin.get();

    pmg::Edge* opp_edge = newOppFEdges[0]->edge;
    Point opp_verts_poses[2];
    opp_verts_poses[0] = opp_edge->verts[0]->pos();
    opp_verts_poses[1] = opp_edge->verts[1]->pos();
    Point main_vert0_pos  = newFFacets[0]->facet->contains(fEdge->edge->verts[0]) ?
                fEdge->edge->verts[0]->pos() : fEdge->edge->verts[1]->pos();
    Point main_vert1_pos  = newFFacets[1]->facet->contains(fEdge->edge->verts[0]) ?
                fEdge->edge->verts[0]->pos() : fEdge->edge->verts[1]->pos();
    Point main_vert0_proj = tva::spatalgs::project(main_vert0_pos, opp_verts_poses[0], opp_verts_poses[1]);
    Point main_vert1_proj = tva::spatalgs::project(main_vert1_pos, opp_verts_poses[0], opp_verts_poses[1]);
    Vec main_vec0 = main_vert0_pos - main_vert0_proj;
    Vec main_vec1 = main_vert1_pos - main_vert1_proj;

    Point adj_opp_pos0 = oppFFacets.first->facet->findVertNot(opp_edge)->pos();
    Point adj_opp_pos1 = oppFFacets.second->facet->findVertNot(opp_edge)->pos();
    Point adj_opp_proj0 = tva::spatalgs::project(adj_opp_pos0, opp_verts_poses[0], opp_verts_poses[1]);
    Point adj_opp_proj1 = tva::spatalgs::project(adj_opp_pos1, opp_verts_poses[0], opp_verts_poses[1]);
    Vec adj_vec0 = adj_opp_pos0 - adj_opp_proj0;
    Vec adj_vec1 = adj_opp_pos1 - adj_opp_proj1;

    double coses[2][2];
    coses[0][0] = Vec::cos(main_vec0, adj_vec0);
    coses[0][1] = Vec::cos(main_vec0, adj_vec1);
    coses[1][0] = Vec::cos(main_vec1, adj_vec0);
    coses[1][1] = Vec::cos(main_vec1, adj_vec1);

    if (coses[0][0] > coses[0][1] && coses[1][1] > coses[1][0])
    {
        newOppFEdges[0]->fillAdjFFacets(newFFacets[0], oppFFacets.first);
        oppFFacets.first->addFEdge(newOppFEdges[0]);

        newOppFEdges[1]->fillAdjFFacets(newFFacets[1], oppFFacets.second);
        oppFFacets.second->addFEdge(newOppFEdges[1]);
        return;
    }

    if (coses[0][1] >= coses[0][0] && coses[1][0] >= coses[1][1])
    {
        newOppFEdges[0]->fillAdjFFacets(newFFacets[0], oppFFacets.second);
        oppFFacets.second->addFEdge(newOppFEdges[0]);

        newOppFEdges[1]->fillAdjFFacets(newFFacets[1], oppFFacets.first);
        oppFFacets.first->addFEdge(newOppFEdges[1]);
        return;
    }

    int best_ofi = -1;
    if (coses[0][0] > coses[0][1] && coses[1][0] > coses[1][1])
        best_ofi = 0;

    if (coses[0][1] >= coses[0][0] && coses[1][1] >= coses[1][0])
        best_ofi = 1;

    int worst_ofi = best_ofi == 0 ? 1 : 0;

    FrSuFacet* opp_ffacets[2] = { oppFFacets.first, oppFFacets.second };
    if (coses[0][best_ofi] > coses[1][best_ofi])
    {
        newOppFEdges[0]->fillAdjFFacets(newFFacets[0], opp_ffacets[best_ofi]);
        opp_ffacets[best_ofi]->addFEdge(newOppFEdges[0]);

        newOppFEdges[1]->fillAdjFFacets(newFFacets[1], opp_ffacets[worst_ofi]);
        opp_ffacets[worst_ofi]->addFEdge(newOppFEdges[1]);
    }
    else
    {
        newOppFEdges[0]->fillAdjFFacets(newFFacets[0], opp_ffacets[worst_ofi]);
        opp_ffacets[worst_ofi]->addFEdge(newOppFEdges[0]);

        newOppFEdges[1]->fillAdjFFacets(newFFacets[1], opp_ffacets[best_ofi]);
        opp_ffacets[best_ofi]->addFEdge(newOppFEdges[1]);
    }

    return;
}


void Crystallite::exhaustFrontCollapse(FrSuEdge *fEdge, FrSuEdge *oppFEdge)
{
    auto fedge_adj_facets = fEdge->getAdjFFacets();
    auto opp_ffacets   = oppFEdge->getAdjFFacets();

    m_innerTetrs.push_back(new Tetr(
        fEdge->edge->verts[0],
        fEdge->edge->verts[1],
        oppFEdge->edge->verts[0],
        oppFEdge->edge->verts[1]));

    auto new_tetr = new Tetr(fEdge->edge->verts[0],
                             fEdge->edge->verts[1],
                             oppFEdge->edge->verts[0],
                             oppFEdge->edge->verts[1]);
    m_innerTetrs.push_back(new_tetr);

//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(Polycrystal::FileType::Obj, "debug.obj");


    std::vector<FrSuEdge*> fedges_to_erase;
    fedges_to_erase.reserve(4);
    for (auto& fedge : opp_ffacets.first->fEdges)
        if (fedge->edge->contains(fEdge->edge->verts[0]) || fedge->edge->contains(fEdge->edge->verts[1]))
            fedges_to_erase.push_back(fedge);

    for (auto& fedge : opp_ffacets.second->fEdges)
        if (fedge->edge->contains(fEdge->edge->verts[0]) || fedge->edge->contains(fEdge->edge->verts[1]))
            fedges_to_erase.push_back(fedge);


    removeFromFront(fedge_adj_facets.first);
    removeFromFront(fedge_adj_facets.second);
    removeFromFront(opp_ffacets.first);
    removeFromFront(opp_ffacets.second);
    removeFromFront(fEdge);
    removeFromFront(oppFEdge);
    for (auto& fedge : fedges_to_erase)
        removeFromFront(fedge);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustFrontCollapse";
#endif
}


void Crystallite::exhaustFrontSplit(FrSuEdge* fEdge, FrSuEdge* oppFEdge)
{
    // This volatile helps to avoid computational error.
    volatile auto adj_f_facets = fEdge->getAdjFFacets();

    Vertex* opp_verts[2];
    opp_verts[0] = adj_f_facets.first->facet->findVertNot(fEdge->edge);
    opp_verts[1] = adj_f_facets.second->facet->findVertNot(fEdge->edge);

    auto opp_edge    = oppFEdge->edge;
    auto opp_ffacets = oppFEdge->getAdjFFacets();

    opp_ffacets.first->removeFEdge(oppFEdge);
    opp_ffacets.second->removeFEdge(oppFEdge);
    removeFromFront(oppFEdge);
    FrSuEdge* new_opp_fedges[2];
    new_opp_fedges[0] = addToFront(opp_edge, false);
    new_opp_fedges[1] = addToFront(opp_edge, false);

    FrSuEdge* new_tetr_fedges[3];
    new_tetr_fedges[0] = nullptr;
    new_tetr_fedges[1] = nullptr;
    new_tetr_fedges[2] = new_opp_fedges[0];

    FrSuFacet* new_ffacets[2];

    for (auto& fedge : adj_f_facets.first->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[0]) &&
            (fedge->edge->contains(opp_edge->verts[0]) ||
             fedge->edge->contains(opp_edge->verts[1])))
        {
            new_tetr_fedges[0] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFacet(adj_f_facets.first);
            fedge->removeAdjFFacet(adj_f_facets.second);
            break;
        }
    }
    for (auto& fedge : adj_f_facets.second->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[0]) &&
            (fedge->edge->contains(opp_edge->verts[0]) ||
             fedge->edge->contains(opp_edge->verts[1])))
        {
            new_tetr_fedges[1] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFacet(adj_f_facets.first);
            fedge->removeAdjFFacet(adj_f_facets.second);
            break;
        }
    }
    new_ffacets[0] = addToFront(new Facet(new_tetr_fedges[0]->edge,
                                          new_tetr_fedges[1]->edge,
                                          new_tetr_fedges[2]->edge));
    new_ffacets[0]->addFEdge(new_tetr_fedges[0]);
    new_ffacets[0]->addFEdge(new_tetr_fedges[1]);
    new_ffacets[0]->addFEdge(new_tetr_fedges[2]);
    new_ffacets[0]->normal = computeNormalInTetr(new_ffacets[0], fEdge->edge);


    new_tetr_fedges[2] = new_opp_fedges[1];

    for (auto& fedge : adj_f_facets.first->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[1]) &&
            (fedge->edge->contains(opp_edge->verts[0]) ||
             fedge->edge->contains(opp_edge->verts[1])))
        {
            new_tetr_fedges[0] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFacet(adj_f_facets.first);
            fedge->removeAdjFFacet(adj_f_facets.second);
            break;
        }
    }
    for (auto& fedge : adj_f_facets.second->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[1]) &&
            (fedge->edge->contains(opp_edge->verts[0]) ||
             fedge->edge->contains(opp_edge->verts[1])))
        {
            new_tetr_fedges[1] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFacet(adj_f_facets.first);
            fedge->removeAdjFFacet(adj_f_facets.second);
            break;
        }
    }
    new_ffacets[1] = addToFront(new Facet(new_tetr_fedges[0]->edge,
                                          new_tetr_fedges[1]->edge,
                                          new_tetr_fedges[2]->edge));
    new_ffacets[1]->addFEdge(new_tetr_fedges[0]);
    new_ffacets[1]->addFEdge(new_tetr_fedges[1]);
    new_ffacets[1]->addFEdge(new_tetr_fedges[2]);
    new_ffacets[1]->normal = computeNormalInTetr(new_ffacets[1], fEdge->edge);


    setFEdgesInFrontSplit(fEdge, new_opp_fedges, new_ffacets, opp_ffacets);

//    m_innerTetrs.push_back(new Tetr(
//        fEdge->edge->verts[0],
//        fEdge->edge->verts[1],
//        opp_verts[0],
//        opp_verts[1]));

    auto new_tetr = new Tetr(fEdge->edge->verts[0],
                             fEdge->edge->verts[1],
                             opp_verts[0],
                             opp_verts[1]);
    m_innerTetrs.push_back(new_tetr);

//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(Polycrystal::FileType::Obj, "debug.obj");


    removeFromFront(adj_f_facets.first);
    removeFromFront(adj_f_facets.second);
    removeFromFront(fEdge);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustFrontSplit";
#endif
}


void Crystallite::exhaustWithoutNewVertOppEdgeExists(FrSuEdge* fEdge, FrSuEdge* oppFEdge)
{
    auto opp_ffacets = oppFEdge->getAdjFFacets();

    FrSuFacet* main_ffacets[3];
    Vertex* main_vert;
    auto fedge_adj_facets = fEdge->getAdjFFacets();
    main_ffacets[0] = std::get<0>(fedge_adj_facets);
    main_ffacets[1] = std::get<1>(fedge_adj_facets);

    if (std::get<0>(opp_ffacets)->facet->contains(fEdge->edge->verts[0]))
    {
        main_ffacets[2] = std::get<0>(opp_ffacets);
        main_vert = fEdge->edge->verts[0];
    }
    else if (std::get<0>(opp_ffacets)->facet->contains(fEdge->edge->verts[1]))
    {
        main_ffacets[2] = std::get<0>(opp_ffacets);
        main_vert = fEdge->edge->verts[1];
    }
    else if (std::get<1>(opp_ffacets)->facet->contains(fEdge->edge->verts[0]))
    {
        main_ffacets[2] = std::get<1>(opp_ffacets);
        main_vert = fEdge->edge->verts[0];
    }
    else
    {
        main_ffacets[2] = std::get<1>(opp_ffacets);
        main_vert = fEdge->edge->verts[1];
    }

    FrSuEdge* new_tetr_fedges[3];
    for (int i = 0; i < 3; i++)
    {
        new_tetr_fedges[i] = main_ffacets[i]->findFEdgeNot(main_vert);
        new_tetr_fedges[i]->refreshAngleData();
        new_tetr_fedges[i]->removeAdjFFacet(main_ffacets[i]);
    }

    auto new_ffacet = addToFront(new Facet(new_tetr_fedges[0]->edge,
                                           new_tetr_fedges[1]->edge,
                                           new_tetr_fedges[2]->edge));
    new_ffacet->addFEdge(new_tetr_fedges[0]);
    new_ffacet->addFEdge(new_tetr_fedges[1]);
    new_ffacet->addFEdge(new_tetr_fedges[2]);
    new_ffacet->normal = computeNormalInTetr(new_ffacet, main_vert->pos());

//    m_innerTetrs.push_back(new Tetr(
//        fEdge->edge->verts[0],
//        fEdge->edge->verts[1],
//        oppFEdge->edge->verts[0],
//        oppFEdge->edge->verts[1]));

    auto new_tetr = new Tetr(fEdge->edge->verts[0],
                             fEdge->edge->verts[1],
                             oppFEdge->edge->verts[0],
                             oppFEdge->edge->verts[1]);
    m_innerTetrs.push_back(new_tetr);

//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(Polycrystal::FileType::Obj, "debug.obj");

    std::vector<FrSuEdge*> erased_fedges;
    erased_fedges.reserve(3);
    for (auto& f_facet : main_ffacets)
    {
        for (auto& fedge : f_facet->fEdges)
        {
            if (std::find(erased_fedges.begin(), erased_fedges.end(), fedge) == erased_fedges.end() &&
                fedge->edge->contains(main_vert))
            {
                removeFromFront(fedge);
                erased_fedges.push_back(fedge);
            }
        }
    }

    for (auto& f_facet : main_ffacets)
        removeFromFront(f_facet);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustWithoutNewVertOppEdgeExists";
#endif
}


void Crystallite::exhaustWithoutNewVertOppEdgeDontExists(FrSuEdge* fEdge)
{
    // This volatile helps to avoid computational error.
    volatile auto adj_f_facets = fEdge->getAdjFFacets();

    Vertex* opp_verts[2];
    opp_verts[0] = adj_f_facets.first->facet->findVertNot(fEdge->edge);
    opp_verts[1] = adj_f_facets.second->facet->findVertNot(fEdge->edge);
    
    FrSuEdge* opp_fedge = addToFront(new Edge(opp_verts[0], opp_verts[1]));

    FrSuEdge* new_tetr_fedges[3];
    new_tetr_fedges[0] = nullptr;
    new_tetr_fedges[1] = nullptr;
    new_tetr_fedges[2] = opp_fedge;

    for (auto& fedge : adj_f_facets.first->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[0]) &&
            (fedge->edge->contains(opp_fedge->edge->verts[0]) ||
             fedge->edge->contains(opp_fedge->edge->verts[1])))
        {
            new_tetr_fedges[0] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFacet(adj_f_facets.first);
//            fedge->removeAdjFFacet(adj_f_facets.second);
            break;
        }
    }
    for (auto& fedge : adj_f_facets.second->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[0]) &&
            (fedge->edge->contains(opp_fedge->edge->verts[0]) ||
             fedge->edge->contains(opp_fedge->edge->verts[1])))
        {
            new_tetr_fedges[1] = fedge;
            fedge->refreshAngleData();
//            fedge->removeAdjFFacet(adj_f_facets.first);
            fedge->removeAdjFFacet(adj_f_facets.second);
            break;
        }
    }
    auto new_ffacet = addToFront(new Facet(new_tetr_fedges[0]->edge,
                                           new_tetr_fedges[1]->edge,
                                           new_tetr_fedges[2]->edge));
    new_ffacet->addFEdge(new_tetr_fedges[0]);
    new_ffacet->addFEdge(new_tetr_fedges[1]);
    new_ffacet->addFEdge(new_tetr_fedges[2]);
    new_ffacet->normal = computeNormalInTetr(new_ffacet, fEdge->edge);

    for (auto& fedge : adj_f_facets.first->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[1]) &&
            (fedge->edge->contains(opp_fedge->edge->verts[0]) ||
             fedge->edge->contains(opp_fedge->edge->verts[1])))
        {
            new_tetr_fedges[0] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFacet(adj_f_facets.first);
            fedge->removeAdjFFacet(adj_f_facets.second);
            break;
        }
    }
    for (auto& fedge : adj_f_facets.second->fEdges)
    {
        if (fedge->edge->contains(fEdge->edge->verts[1]) &&
            (fedge->edge->contains(opp_fedge->edge->verts[0]) ||
             fedge->edge->contains(opp_fedge->edge->verts[1])))
        {
            new_tetr_fedges[1] = fedge;
            fedge->refreshAngleData();
            fedge->removeAdjFFacet(adj_f_facets.first);
            fedge->removeAdjFFacet(adj_f_facets.second);
            break;
        }
    }
    new_ffacet = addToFront(new Facet(new_tetr_fedges[0]->edge,
                                      new_tetr_fedges[1]->edge,
                                      new_tetr_fedges[2]->edge));
    new_ffacet->addFEdge(new_tetr_fedges[0]);
    new_ffacet->addFEdge(new_tetr_fedges[1]);
    new_ffacet->addFEdge(new_tetr_fedges[2]);
    new_ffacet->normal = computeNormalInTetr(new_ffacet, fEdge->edge);

//    m_innerTetrs.push_back(new Tetr(
//        fEdge->edge->verts[0],
//        fEdge->edge->verts[1],
//        opp_verts[0],
//        opp_verts[1]));

    auto new_tetr = new Tetr(fEdge->edge->verts[0],
                             fEdge->edge->verts[1],
                             opp_verts[0],
                             opp_verts[1]);
    m_innerTetrs.push_back(new_tetr);

//    std::cout << std::endl << fEdge->edge->magnitude() << ' ' << fEdge->computeAngle() * 180.0 / M_PI;
//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(Polycrystal::FileType::Obj, "debug.obj");

    
    removeFromFront(adj_f_facets.first);
    removeFromFront(adj_f_facets.second);
    removeFromFront(fEdge);

#ifdef DEV_DEBUG
    std::cout << std::endl << "exhaustWithoutNewVertOppEdgeDontExists";
#endif
}


void Crystallite::exhaustWithoutNewVert(FrSuEdge* fEdge, bool oppEdgeExistence, FrSuEdge* oppFEdge)
{
    FrSuEdge* opp_fedge = nullptr;
    if (oppEdgeExistence && oppFEdge)
        opp_fedge = oppFEdge;
    else if (oppEdgeExistence)
        opp_fedge = fEdge->findOppEdge();

    if (/*(oppFEdge && oppEdgeExistence) ||*/
        opp_fedge)
    {
        if (frontCollapseCheck(fEdge, opp_fedge))
        {
            exhaustFrontCollapse(fEdge, opp_fedge);
        }
        else if (frontSplitCheck(fEdge, opp_fedge))
        {
            exhaustFrontSplit(fEdge, opp_fedge);
        }
        else
        {
            exhaustWithoutNewVertOppEdgeExists(fEdge, opp_fedge);
        }
    }
    else
    {
        exhaustWithoutNewVertOppEdgeDontExists(fEdge);
    }
}




bool Crystallite::tryComputeNewVertPosType3(FrSuFacet* fFacet, Vec& out_pos)
{
    FrSuEdge* main_fedges[2];
    main_fedges[0] = fFacet->fEdges[0];
    main_fedges[1] = fFacet->fEdges[1];
    Edge* main_edges[2];
    main_edges[0] = main_fedges[0]->edge;
    main_edges[1] = main_fedges[1]->edge;

    auto main_vert = main_edges[0]->verts[0] == main_edges[1]->verts[0] ?
        main_edges[0]->verts[0] :
        (main_edges[0]->verts[0] == main_edges[1]->verts[1] ?
            main_edges[0]->verts[0] :
            main_edges[0]->verts[1]);
    auto third_edge = fFacet->facet->findEdgeNot(main_vert);
    auto third_f_edge = fFacet->findFEdge(third_edge);

    auto v0 = main_edges[0]->verts[0] == main_vert ? main_edges[0]->verts[1] : main_edges[0]->verts[0];
    auto v1 = main_edges[1]->verts[0] == main_vert ? main_edges[1]->verts[1] : main_edges[1]->verts[0];
    auto v2 = main_vert;

    auto adj_f_facets = main_fedges[0]->getAdjFFacets();
    auto fn0 = std::get<0>(adj_f_facets) == fFacet ? std::get<1>(adj_f_facets) : std::get<0>(adj_f_facets);
    adj_f_facets = main_fedges[1]->getAdjFFacets();
    auto fn1 = std::get<0>(adj_f_facets) == fFacet ? std::get<1>(adj_f_facets) : std::get<0>(adj_f_facets);
    adj_f_facets = third_f_edge->getAdjFFacets();
    auto fn2 = std::get<0>(adj_f_facets) == fFacet ? std::get<1>(adj_f_facets) : std::get<0>(adj_f_facets);

    const Vec& v0_pos = v0->pos();
    const Vec& v1_pos = v1->pos();
    const Vec& v2_pos = v2->pos();

    Vec n_m = fFacet->normal;
    Vec n_n0 = fn0->normal;
    Vec n_n1 = fn1->normal;
    Vec n_n2 = fn2->normal;

    Vec e_mn0 = (n_m + n_n0).normalize();
    Vec e_mn1 = (n_m + n_n1).normalize();
    Vec e_mn2 = n_m + n_n2;

    Vec np_mn0 = Vec::cross(v2_pos - v0_pos, e_mn0);
    Vec np_mn1 = Vec::cross(v2_pos - v1_pos, e_mn1);

    Vec e = Vec::cross(np_mn0, np_mn1);

    Vec new_pos = tva::spatalgs::lineIntersectPlane(v2_pos, e, v0_pos, v1_pos, v0_pos + e_mn2);

    Vec v0_to_np = new_pos - v0_pos;
    Vec v1_to_np = new_pos - v1_pos;
    Vec v2_to_np = new_pos - v2_pos;
    double sum_magn0 =  fn0->facet->edges[0]->magnitude()
                      + fn0->facet->edges[1]->magnitude()
                      + fn0->facet->edges[2]->magnitude();
    double sum_magn1 =  fn1->facet->edges[0]->magnitude()
                      + fn1->facet->edges[1]->magnitude()
                      + fn1->facet->edges[2]->magnitude();
    double sum_magn2 =  fn2->facet->edges[0]->magnitude()
                      + fn2->facet->edges[1]->magnitude()
                      + fn2->facet->edges[2]->magnitude();
    double sum_magn3 =  fFacet->facet->edges[0]->magnitude()
                      + fFacet->facet->edges[1]->magnitude()
                      + fFacet->facet->edges[2]->magnitude();
    double av_magn = (sum_magn0 + sum_magn1 + sum_magn2 + sum_magn3) / 12.0;
    if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
//        doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * av_magn) ||
        !vertInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        return false;
    }

//    std::cout << std::endl << "Type3";

    out_pos = new_pos;
    return true;
}


bool Crystallite::tryComputeNewVertPosType2(FrSuFacet* fFacet, Vec& out_pos, int smallAngleIndex0, int smallAngleIndex1)
{
    FrSuEdge* main_fedges[2];
    main_fedges[0] = fFacet->fEdges[smallAngleIndex0];
    main_fedges[1] = fFacet->fEdges[smallAngleIndex1];
    Edge* main_edges[2];
    main_edges[0] = main_fedges[0]->edge;
    main_edges[1] = main_fedges[1]->edge;
    auto main_vert = main_edges[0]->verts[0] == main_edges[1]->verts[0] ?
        main_edges[0]->verts[0] :
        (main_edges[0]->verts[0] == main_edges[1]->verts[1] ?
            main_edges[0]->verts[0] :
            main_edges[0]->verts[1]);

    auto v0 = main_edges[0]->verts[0] == main_vert ? main_edges[0]->verts[1] : main_edges[0]->verts[0];
    auto v1 = main_edges[1]->verts[0] == main_vert ? main_edges[1]->verts[1] : main_edges[1]->verts[0];
    auto v2 = main_vert;

    auto adj_f_facets = main_fedges[0]->getAdjFFacets();
    auto fn0 = std::get<0>(adj_f_facets) == fFacet ? std::get<1>(adj_f_facets) : std::get<0>(adj_f_facets);
    auto vn0 = fn0->facet->findVertNot(main_fedges[0]->edge);
    adj_f_facets = main_fedges[1]->getAdjFFacets();
    auto fn1 = std::get<0>(adj_f_facets) == fFacet ? std::get<1>(adj_f_facets) : std::get<0>(adj_f_facets);
    auto vn1 = fn1->facet->findVertNot(main_fedges[1]->edge);

    Point v0_pos = v0->pos();
    Point v1_pos = v1->pos();
    Point v2_pos = v2->pos();
    Point vn0_pos = vn0->pos();
    Point vn1_pos = vn1->pos();

    Vec n_m = fFacet->normal;
    Vec n_n0 = fn0->normal;
    Vec n_n1 = fn1->normal;
    
    Vec e_mn0 = (n_m + n_n0).normalize();
    Vec e_mn1 = (n_m + n_n1).normalize();

    Vec np_mn0 = Vec::cross(v2_pos - v0_pos, e_mn0);
    Vec np_mn1 = Vec::cross(v2_pos - v1_pos, e_mn1);

    Vec e = Vec::cross(np_mn0, np_mn1).normalize();
    if (Vec::dot(e, n_m) < 0.0) e *= -1.0;

//    double f_facet_area = fFacet->facet->computeArea();
//    double sp = SQRT3_2 * 0.5 * m_preferredLength * m_preferredLength;
//    double sf0 = 0.5 * (f_facet_area + fn0->facet->computeArea());
//    double sf1 = 0.5 * (f_facet_area + fn1->facet->computeArea());
//    double raw_deform0 = K_D * (sp - sf0);
//    double raw_deform1 = K_D * (sp - sf1);
//    double deform0 = raw_deform0 < sf0 * K_MAXD ? raw_deform0 : sf0 * K_MAXD;
//    double deform1 = raw_deform1 < sf1 * K_MAXD ? raw_deform1 : sf1 * K_MAXD;
//    double sc0 = sf0 + deform0;
//    double sc1 = sf1 + deform1;
//    double x0_2 = sc0 / Vec::cross(v2->pos() - v0->pos(), e).magnitude();
//    double x1_2 = sc1 / Vec::cross(v2->pos() - v1->pos(), e).magnitude();
//    Vec new_pos = v2_pos + (x0_2 + x1_2) * e;
    double l0 = (v1_pos - v0_pos).magnitude();
    double l1 = (vn0_pos - v0_pos).magnitude();
    double l2 = (vn0_pos - v2_pos).magnitude();
    double l3 = (vn1_pos - v1_pos).magnitude();
    double l4 = (vn1_pos - v2_pos).magnitude();
    double av_magn = (main_edges[0]->magnitude() + main_edges[1]->magnitude() + l0 + l1 + l2 + l3 + l4) / 7.0;
    double raw_deform = K_D * (m_preferredLength - av_magn);
    double deform = raw_deform < av_magn * K_MAXD ? raw_deform : av_magn * K_MAXD;
    double magn_d = av_magn + deform;
    Vec new_pos = v2_pos + magn_d * e;

    Vec v0_to_np = new_pos - v0_pos;
    Vec v1_to_np = new_pos - v1_pos;
    Vec v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
//        doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * magn_d) ||
        !vertInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
//        x0_2 = sf0 / Vec::cross(v2->pos() - v0->pos(), e).magnitude();
//        x1_2 = sf1 / Vec::cross(v2->pos() - v1->pos(), e).magnitude();
        new_pos = new_pos = v2_pos + av_magn * e;
        v0_to_np = new_pos - v0_pos;
        v1_to_np = new_pos - v1_pos;
        v2_to_np = new_pos - v2_pos;
        if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
            segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
            segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
//            doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * av_magn) ||
            !vertInsideFrontCheck(new_pos) ||
            facetIntersectionCheck(v0, v1, new_pos) ||
            facetIntersectionCheck(v0, v2, new_pos) ||
            facetIntersectionCheck(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type2";

    out_pos = new_pos;
    return true;
}


bool Crystallite::tryComputeNewVertPosType1(FrSuFacet* fFacet, Vec& out_pos, int smallAngleIndex)
{
    auto main_f_edge = fFacet->fEdges[smallAngleIndex];
    auto main_edge = fFacet->fEdges[smallAngleIndex]->edge;
    Vec main_edge_poses[2];
    main_edge_poses[0] = main_edge->verts[0]->pos();
    main_edge_poses[1] = main_edge->verts[1]->pos();

    auto v0 = main_edge->verts[0];
    auto v1 = main_edge->verts[1];
    auto v2 = fFacet->facet->findVertNot(main_edge);

    auto adj_f_facets = main_f_edge->getAdjFFacets();
    auto fn = std::get<0>(adj_f_facets) == fFacet ? std::get<1>(adj_f_facets) : std::get<0>(adj_f_facets);
    auto vn = fn->facet->findVertNot(main_edge);

    Point v0_pos = v0->pos();
    Point v1_pos = v1->pos();
    Point v2_pos = v2->pos();
    Point vn_pos = vn->pos();

    Point v2pr = tva::spatalgs::project(v2_pos, v0_pos, v1_pos);
    Point vnpr = tva::spatalgs::project(vn_pos, v0_pos, v1_pos);

    Point c = 0.25 * (v0_pos + v1_pos + v2pr + vnpr);
    Vec   e = (fFacet->normal + fn->normal).normalize();
    double me_magn = main_edge->magnitude();
    double l0 = (v2_pos - v0_pos).magnitude();
    double l1 = (v2_pos - v1_pos).magnitude();
    double l2 = (vn_pos - v0_pos).magnitude();
    double l3 = (vn_pos - v1_pos).magnitude();
    double av_magn = 0.2 * (me_magn + l0 + l1 + l2 + l3);
    double v0_c_dist = (c - v0_pos).magnitude();
    double raw_deform = K_D * (m_preferredLength - av_magn);
    double deform = raw_deform < av_magn * K_MAXD ? raw_deform : av_magn * K_MAXD;
    double magn_d = av_magn + deform;
    Point new_pos = c + sqrt(magn_d * magn_d - v0_c_dist * v0_c_dist) * e;

    Vec v0_to_np = new_pos - v0_pos;
    Vec v1_to_np = new_pos - v1_pos;
    Vec v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
//        doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * magn_d) ||
        !vertInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        new_pos = c + sqrt(av_magn * av_magn - v0_c_dist * v0_c_dist) * e;
        v0_to_np = new_pos - v0_pos;
        v1_to_np = new_pos - v1_pos;
        v2_to_np = new_pos - v2_pos;
        if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
            segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
            segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
//            doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * av_magn) ||
            !vertInsideFrontCheck(new_pos) ||
            facetIntersectionCheck(v0, v1, new_pos) ||
            facetIntersectionCheck(v0, v2, new_pos) ||
            facetIntersectionCheck(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type1";

    out_pos = new_pos;
    return true;
}


bool Crystallite::tryComputeNewVertPosType0(FrSuFacet* fFacet, Vec& out_pos)
{
    double av_magn = 0.3333333333333 * (
          fFacet->facet->edges[0]->magnitude()
        + fFacet->facet->edges[1]->magnitude()
        + fFacet->facet->edges[2]->magnitude());
    double raw_deform = K_D * (m_preferredLength - av_magn);
    double deform = raw_deform < av_magn * K_MAXD ? raw_deform : av_magn * K_MAXD;
    double magn_d = av_magn + deform;
    Vec new_pos = fFacet->computeCenter() + sqrt(magn_d * magn_d - 0.33333333333 * av_magn * av_magn) * fFacet->normal;

    auto v0 = fFacet->facet->edges[0]->verts[0];
    auto v1 = fFacet->facet->edges[0]->verts[1];
    auto v2 = fFacet->facet->findVertNot(fFacet->facet->edges[0]);
    const Vec& v0_pos = v0->pos();
    const Vec& v1_pos = v1->pos();
    const Vec& v2_pos = v2->pos();
    Vec v0_to_np = new_pos - v0_pos;
    Vec v1_to_np = new_pos - v1_pos;
    Vec v2_to_np = new_pos - v2_pos;
    if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
        segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
        segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
//        doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * magn_d) ||
        !vertInsideFrontCheck(new_pos) ||
        facetIntersectionCheck(v0, v1, new_pos) ||
        facetIntersectionCheck(v0, v2, new_pos) ||
        facetIntersectionCheck(v1, v2, new_pos))
    {
        new_pos = fFacet->computeCenter() + sqrt(av_magn * av_magn - 0.33333333333 * av_magn * av_magn) * fFacet->normal;
        v0_to_np = new_pos - v0_pos;
        v1_to_np = new_pos - v1_pos;
        v2_to_np = new_pos - v2_pos;
        if (segmentFrontIntersectionCheck(v0_pos + FROM_VERT_COEF * v0_to_np, new_pos + NOT_TOO_CLOSE * v0_to_np) ||
            segmentFrontIntersectionCheck(v1_pos + FROM_VERT_COEF * v1_to_np, new_pos + NOT_TOO_CLOSE * v1_to_np) ||
            segmentFrontIntersectionCheck(v2_pos + FROM_VERT_COEF * v2_to_np, new_pos + NOT_TOO_CLOSE * v2_to_np) ||
//            doesFrontIntersectSphere(new_pos, NOT_TOO_CLOSE * av_magn) ||
            !vertInsideFrontCheck(new_pos) ||
            facetIntersectionCheck(v0, v1, new_pos) ||
            facetIntersectionCheck(v0, v2, new_pos) ||
            facetIntersectionCheck(v1, v2, new_pos))
            return false;
    }

//    std::cout << std::endl << "Type0";

    out_pos = new_pos;
    return true;
}


bool Crystallite::tryComputeNewVertPos(FrSuFacet* fFacet, Vec& out_pos)
{
    double angs_coses[3]
    { 
        fFacet->fEdges[0]->angleExCos(),
        fFacet->fEdges[1]->angleExCos(),
        fFacet->fEdges[2]->angleExCos()
    };
    int indexes[3];
    int small_angs_num = 0;
    if (angs_coses[0] > cosDeg<140>) indexes[small_angs_num++] = 0;
    if (angs_coses[1] > cosDeg<140>) indexes[small_angs_num++] = 1;
    if (angs_coses[2] > cosDeg<140>) indexes[small_angs_num++] = 2;

    switch (small_angs_num)
    {
    case 0: return tryComputeNewVertPosType0(fFacet, out_pos);
    case 1: return tryComputeNewVertPosType1(fFacet, out_pos, indexes[0]);
    case 2: return tryComputeNewVertPosType2(fFacet, out_pos, indexes[0], indexes[1]);
    case 3: return tryComputeNewVertPosType3(fFacet, out_pos);
    }

    return true;
}




double Crystallite::sqr4FacetArea(const FrSuFacet* fFacet) const
{
    Vec vec0 = *fFacet->facet->edges[0]->verts[1] - *fFacet->facet->edges[0]->verts[0];
    Vec vec1 = *fFacet->facet->edges[1]->verts[1] - *fFacet->facet->edges[1]->verts[0];
    return Vec::cross(vec0, vec1).sqrMagnitude();
}


FrSuFacet* Crystallite::chooseFacetForExhaustionWithNewVert(FrSuEdge* fEdge)
{
    auto adj_f_facets = fEdge->getAdjFFacets();

    return sqr4FacetArea(std::get<0>(adj_f_facets)) < sqr4FacetArea(std::get<1>(adj_f_facets)) ?
//    return std::get<0>(adj_f_facets)->computeQuality() < std::get<1>(adj_f_facets)->computeQuality() ?
        std::get<0>(adj_f_facets) :
        std::get<1>(adj_f_facets);
}


void Crystallite::exhaustWithNewVert(FrSuFacet* fFacet, const Vec& vertPos)
{
    Vertex* new_vert = new Vertex(vertPos);
    m_innerVerts.push_back(new_vert);

    FrSuEdge* new_tetr_fedges[6];
    new_tetr_fedges[0] = fFacet->fEdges[0];
    auto far_vert = fFacet->facet->findVertNot(new_tetr_fedges[0]->edge);
    new_tetr_fedges[1] = fFacet->findFEdge(new_tetr_fedges[0]->edge->verts[1], far_vert);
    new_tetr_fedges[2] = fFacet->findFEdge(new_tetr_fedges[0]->edge->verts[0], far_vert);
    new_tetr_fedges[3] = addToFront(new Edge(new_tetr_fedges[0]->edge->verts[0], new_vert));
    new_tetr_fedges[4] = addToFront(new Edge(new_tetr_fedges[0]->edge->verts[1], new_vert));
    new_tetr_fedges[5] = addToFront(new Edge(far_vert, new_vert));

    new_tetr_fedges[0]->refreshAngleData();
    new_tetr_fedges[1]->refreshAngleData();
    new_tetr_fedges[2]->refreshAngleData();

    new_tetr_fedges[0]->removeAdjFFacet(fFacet);
    new_tetr_fedges[1]->removeAdjFFacet(fFacet);
    new_tetr_fedges[2]->removeAdjFFacet(fFacet);
    
    auto new_ffacet = addToFront(new Facet(new_tetr_fedges[0]->edge,
                                           new_tetr_fedges[3]->edge,
                                           new_tetr_fedges[4]->edge));
    new_ffacet->addFEdge(new_tetr_fedges[0]);
    new_ffacet->addFEdge(new_tetr_fedges[3]);
    new_ffacet->addFEdge(new_tetr_fedges[4]);
    new_ffacet->normal = computeNormalInTetr(new_ffacet, far_vert->pos());

    new_ffacet = addToFront(new Facet(new_tetr_fedges[2]->edge,
                                      new_tetr_fedges[3]->edge,
                                      new_tetr_fedges[5]->edge));
    new_ffacet->addFEdge(new_tetr_fedges[2]);
    new_ffacet->addFEdge(new_tetr_fedges[3]);
    new_ffacet->addFEdge(new_tetr_fedges[5]);
    new_ffacet->normal = computeNormalInTetr(new_ffacet, new_tetr_fedges[0]->edge->verts[1]->pos());

    new_ffacet = addToFront(new Facet(new_tetr_fedges[1]->edge,
                                       new_tetr_fedges[4]->edge,
                                       new_tetr_fedges[5]->edge));
    new_ffacet->addFEdge(new_tetr_fedges[1]);
    new_ffacet->addFEdge(new_tetr_fedges[4]);
    new_ffacet->addFEdge(new_tetr_fedges[5]);
    new_ffacet->normal = computeNormalInTetr(new_ffacet, new_tetr_fedges[0]->edge->verts[0]->pos());

//    m_innerTetrs.push_back(new Tetr(
//        new_tetr_fedges[0]->edge->verts[0],
//        new_tetr_fedges[0]->edge->verts[1],
//        new_tetr_fedges[5]->edge->verts[0],
//        new_tetr_fedges[5]->edge->verts[1]));

    auto new_tetr = new Tetr(new_tetr_fedges[0]->edge->verts[0],
                             new_tetr_fedges[0]->edge->verts[1],
                             new_tetr_fedges[5]->edge->verts[0],
                             new_tetr_fedges[5]->edge->verts[1]);
    m_innerTetrs.push_back(new_tetr);
//    if (new_tetr->computeQuality() < 1e-2)
//        m_polycr->output(Polycrystal::FileType::Obj, "debug.obj");

    removeFromFront(fFacet);
}




bool Crystallite::tryExhaustWithoutNewVert(FrSuEdge* fEdge, bool oppEdgeExistence, FrSuEdge* oppEdge)
{
    if (parallelFacetsCheck(fEdge) ||
        edgeIntersectionCheck(fEdge) ||
        facetsIntersectionCheck(fEdge) ||
        anyVertInsidePotentialTetrCheck(fEdge))
        return false;

    exhaustWithoutNewVert(fEdge, oppEdgeExistence, oppEdge);
    return true;
}


bool Crystallite::tryExhaustWithNewVert(FrSuEdge* frontEdge)
{
    if (parallelFacetsCheck(frontEdge))
        return false;

    auto exhaust_f_facet = chooseFacetForExhaustionWithNewVert(frontEdge);
    Vec new_vert_pos;
    if (!tryComputeNewVertPos(exhaust_f_facet, new_vert_pos))
        return false;

    exhaustWithNewVert(exhaust_f_facet, new_vert_pos);
    return true;
}




bool Crystallite::isFrontExhausted()
{
    if (m_frontFacets.size() == 0 &&
        m_frontEdges.size() == 0)
        return true;

    if ((m_frontFacets.size() == 0 && m_frontEdges.size() > 0) ||
        (m_frontFacets.size() > 0 && m_frontEdges.size() == 0))
        throw std::logic_error("Error in Crystallite::isFrontExhausted. Front wasn't correctly exhausted.");

    return false;
}


void Crystallite::processAngles()
{
#ifdef DEV_DEBUG
    int debug_i = 0;
#endif
    double max_compl = std::numeric_limits<double>::max();
    for (FrSuEdge* cur_fedge = currentFrontEdge(max_compl);; cur_fedge = currentFrontEdge(max_compl))
    {
        if (!cur_fedge)
            throw std::logic_error("pmg::Crystallite::currentFrontEdge returned nullptr");
        
        if (exhaustWithoutNewVertPriorityPredicate(cur_fedge))
        {
            if (!tryExhaustWithoutNewVert(cur_fedge))
            {
                max_compl = cur_fedge->complexity();
                continue;
            }
        }
        else if (exhaustWithNewVertPriorityPredicate(cur_fedge))
        {
            if (!tryExhaustWithNewVert(cur_fedge))
            {
                max_compl = cur_fedge->complexity();
                continue;
            }
        }
        else
        {
            FrSuFacet* exhaust_from_f_facet = nullptr;
            Vec* new_vert_pos = nullptr;
            switch (computeExhaustionTypeQualityPriority(cur_fedge, exhaust_from_f_facet, new_vert_pos))
            {
            case ExhaustType::WithoutNewVert:
                exhaustWithoutNewVert(cur_fedge);
                break;

            case ExhaustType::WithNewVert:
                if (new_vert_pos)
                {
                    exhaustWithNewVert(exhaust_from_f_facet, *new_vert_pos);
                    delete new_vert_pos;
                }
                else
                {
                    if (!tryExhaustWithNewVert(cur_fedge))
                    {
                        max_compl = cur_fedge->complexity();
                        continue;
                    }
                }
                break;

            case ExhaustType::DontExhaust:
                max_compl = cur_fedge->complexity();
#ifdef DEV_DEBUG
                std::cout << std::endl << "DontExhaust";
#endif
                continue;
            }
        }
        max_compl = std::numeric_limits<double>::max();

#ifdef DEV_DEBUG
//        if (debug_i++ >= 1200)
//        m_polycr->output(Polycrystal::FileType::Obj, "debug.obj");
//        if (debug_i++ >= 0)
//            std::cout << std::endl << debug_i - 1;

//        debug();
#endif

        if (isFrontExhausted())
            return;
    }
}


void Crystallite::debug()
{
    Vec accum;
    for (auto& vert : m_innerVerts)
        accum += vert->pos();

    for (auto& svert : m_shellVerts)
        accum += svert->attachedVert->pos();

    for (auto& sedge : m_shellEdges)
        for (auto& vert : sedge->innerVerts())
            accum += vert->pos();

    for (auto& sfacet : m_shellFacets)
        for (auto& vert : sfacet->innerVerts())
            accum += vert->pos();

    std::cout << "\n{ " << accum.coors[0] + accum.coors[1] + accum.coors[2] << " }";
}




void Crystallite::generateMesh(double preferredLen)
{
    m_preferredLength = preferredLen;
    initializeFront();
    computeFrontNormals();
    processAngles();
//    if (globalIntersectionCheck())
//        throw std::logic_error("Intersection error.\npmg::Crystallite3::globalIntersectionCheck returned true.");
    smoothMesh(20);
}


bool Crystallite::globalIntersectionCheck()
{
    for (auto& edge : m_innerEdges)
        if (edgeGlobalIntersectionCheck(edge))
            return true;

    return false;
}




void Crystallite::smoothMesh(unsigned iterationsNum)
{
    for (unsigned i = 0; i < iterationsNum; i++)
    {
        for (auto& vert : m_innerVerts)
        {
            Vec shift;
            int delta_shifts_num = 0;
            for (auto& edge : m_innerEdges)
            {
                if (vert == edge->verts[0])
                {
                    shift += *edge->verts[1] - *vert;
                    delta_shifts_num++;
                }
                else if (vert == edge->verts[1])
                {
                    shift += *edge->verts[0] - *vert;
                    delta_shifts_num++;
                }
            }
            shift /= delta_shifts_num;
            vert->setPos(vert->pos() + shift);
        }
    }
}


void Crystallite::smoothNotFinisedMesh(unsigned iterationsNum)
{
    for (unsigned i = 0; i < iterationsNum; i++)
        for (auto &vert : m_innerVerts)
            smoothAroundFrontVert(vert);
}


void Crystallite::smoothFront(unsigned iterationsNum)
{
    for (unsigned i = 0; i < iterationsNum; i++)
        for (auto &vert : m_innerVerts)
            smoothAroundFrontVert(vert);
}


void Crystallite::smoothAroundFrontVert(Vertex* fVert)
{
    if (fVert->belongsToShellFacet ||
        fVert->belongsToShellEdge ||
        fVert->belongsToShellVertex)
        return;

    Vec shift;
    int delta_shifts_num = 0;
    for (auto &edge : m_innerEdges)
    {
        Vec d_shift;
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
    fVert->setPos(fVert->pos() + shift);
}


std::pair<double, double> Crystallite::analyzeMeshQuality()
{
    size_t simps_num = 0;
    double av_q = 0.0;
    double min_q = 1.0;
    for (auto &tetr : m_innerTetrs)
    {
        double q = tetr->computeQuality();
        av_q += q;
        simps_num++;
        if (q < min_q)
            min_q = q;
    }
    av_q /= simps_num;

    return { min_q, av_q };
}




Crystallite::Crystallite() {}


Crystallite::Crystallite(Polycrystal* polycr)
{
    m_polycr = polycr;
}


Crystallite::~Crystallite()
{
    for (auto& tetr : m_innerTetrs)
        delete tetr;
    for (auto& facet : m_innerFacets)
        delete facet;
    for (auto& edge : m_innerEdges)
        delete edge;
    for (auto& vert : m_innerVerts)
        delete vert;
}




double Crystallite::preferredLength() const
{
    return m_preferredLength;
}


void Crystallite::addToShell(const shell::Facet* shellFacet)
{
    m_shellFacets.push_back(const_cast<shell::Facet*>(shellFacet));
}


void Crystallite::addToShell(const shell::Edge* shellEdge)
{
    m_shellEdges.push_back(const_cast<shell::Edge*>(shellEdge));
}


void Crystallite::addToShell(const shell::Vertex* shellVert)
{
    m_shellVerts.push_back(const_cast<shell::Vertex*>(shellVert));
}




bool Crystallite::shellContains(const shell::Facet* shellFacet) const
{
    return std::find(m_shellFacets.begin(), m_shellFacets.end(), shellFacet) != m_shellFacets.end();
}


bool Crystallite::shellContains(const shell::Edge* shellEdge) const
{
    return std::find(m_shellEdges.begin(), m_shellEdges.end(), shellEdge) != m_shellEdges.end();
}


bool Crystallite::shellContains(const shell::Vertex* shellVert) const
{
    return std::find(m_shellVerts.begin(), m_shellVerts.end(), shellVert) != m_shellVerts.end();
}




const std::vector<Tetr*>& Crystallite::innerTetrs() const
{
    return m_innerTetrs;
}


const std::vector<Facet*>& Crystallite::innerFacets() const
{
    return m_innerFacets;
}


const std::vector<Vertex*>& Crystallite::innerVerts() const
{
    return m_innerVerts;
}




const std::list<FrSuFacet*>& Crystallite::frontFacets() const
{
    return m_frontFacets;
}


const std::list<FrSuEdge*>& Crystallite::frontEdges() const
{
    return m_frontEdges;
}
