#include "Facet3.h"
#include "SpatialAlgs.h"



tva::Vec3 Facet3::computeCenter()
{
    return 0.3333333333333333 * (
        edges[0]->verts[0]->getPos() +
        edges[0]->verts[1]->getPos() +
        findVertNotIncludedInEdge(edges[0])->getPos());
}


double Facet3::computeQuality()
{
    return findShortestEdge()->sqrMagnitude() / findLongestEdge()->sqrMagnitude();
}




Edge3* Facet3::intersectAlongAnEdge(const Facet3* facet0, const Facet3* facet1)
{
    int inters = 0;
    Edge3* res = nullptr;
    for (auto& edge0 : facet0->edges)
        for (auto& edge1 : facet1->edges)
            if (edge0 == edge1)
            {
                inters++;
                res = edge0;
                break;
            }

    return inters == 1 ? res : nullptr;
}




bool Facet3::intersectsBy(const tva::Point3& origin, const tva::Vec3& dir)
{
    return tva::spatialalgs::isRayIntersectTriangle(
        origin, dir,
        edges[0]->verts[0]->getPos(), 
        edges[0]->verts[1]->getPos(), 
        findVertNotIncludedInEdge(edges[0])->getPos());
}


Vertex3* Facet3::findVertNotIncludedInEdge(const Edge3* edge) const
{
    for (auto &facet_edge : edges)
    {
        if (edge != facet_edge)
        {
            if (!edge->contains(facet_edge->verts[0]))
                return facet_edge->verts[0];
            else
                return facet_edge->verts[1];
        }
    }

    return nullptr;
}


Edge3* Facet3::findEdgeNotContainingVert(const Vertex3* vert) const
{
    for (auto &edge : edges)
        if (!edge->contains(vert))
            return edge;

    return nullptr;
}


Edge3* Facet3::findEdge(const Vertex3* vert0, const Vertex3* vert1)
{
    for (auto &edge : edges)
        if ((edge->verts[0] == vert0 &&
             edge->verts[1] == vert1) ||
            (edge->verts[0] == vert1 &&
             edge->verts[1] == vert0))
            return edge;
    
    return nullptr;
}


Edge3* Facet3::findShortestEdge()
{
    if (edges[0]->sqrMagnitude() < edges[1]->sqrMagnitude())
    {
        if (edges[0]->sqrMagnitude() < edges[2]->sqrMagnitude())
            return edges[0];
        else
            return edges[2];
    }
    else
    {
        if (edges[1]->sqrMagnitude() < edges[2]->sqrMagnitude())
            return edges[1];
        else
            return edges[2];
    }
}


Edge3* Facet3::findLongestEdge()
{
    if (edges[0]->sqrMagnitude() > edges[1]->sqrMagnitude())
    {
        if (edges[0]->sqrMagnitude() > edges[2]->sqrMagnitude())
            return edges[0];
        else
            return edges[2];
    }
    else
    {
        if (edges[1]->sqrMagnitude() > edges[2]->sqrMagnitude())
            return edges[1];
        else
            return edges[2];
    }
}




bool Facet3::contains(const Edge3* edge) const
{
    for (auto& edge_ : edges)
        if (edge_ == edge)
            return true;

    return false;
}


bool Facet3::contains(const Vertex3* vert) const
{
    for (auto& edge : edges)
        if (edge->contains(vert))
            return true;

    return false;
}




Facet3::Facet3()
{
    edges[0] = nullptr;
    edges[1] = nullptr;
    edges[2] = nullptr;
}


Facet3::Facet3(const Edge3* edge0, const Edge3* edge1, const Edge3* edge2)
{
    edges[0] = (Edge3*)edge0;
    edges[1] = (Edge3*)edge1;
    edges[2] = (Edge3*)edge2;
}


Facet3::Facet3(const Vertex3* vert0, const Vertex3* vert1, const Vertex3* vert2)
{
    edges[0] = new Edge3(vert0, vert1);
    edges[1] = new Edge3(vert1, vert2);
    edges[2] = new Edge3(vert2, vert0);
}


Facet3::~Facet3() {}




tva::Vec3 FrontFacet3::getNormal()
{
    return m_normal;
}


void FrontFacet3::setNormal(const tva::Vec3& vec)
{
    m_normal = vec;
}


tva::Vec3 FrontFacet3::computeNormal()
{
    tva::Vec3 center = computeCenter();
    tva::Vec3 third_pos = facet->findVertNotIncludedInEdge(facet->edges[0])->getPos();
    tva::Vec3 normal = tva::Vec3::crossProduct(
        facet->edges[0]->verts[0]->getPos() - third_pos,
        facet->edges[0]->verts[1]->getPos() - third_pos).normalize();

    // I don't know why it led to better(correct) result.
    tva::Vec3 test_normal_correct_intersect = normal + tva::Vec3(2.1632737147, 1.488313178, -0.71123534278) * 1e-3;
        
    int intersects_num = 0;
    for (auto& f_facet : m_relatedCrys->getFrontFacets())
    {
        if (f_facet == this)
            continue;

        if (f_facet->facet->intersectsBy(center, /*normal*/test_normal_correct_intersect))
            intersects_num++;
    }

    return m_normal = intersects_num % 2 == 1 ? normal : -normal;
}




tva::Vec3 FrontFacet3::computeCenter()
{
    return facet->computeCenter();
}


double FrontFacet3::computeQuality()
{
    return facet->computeQuality();
}




FrontFacet3::FrontFacet3(const Crystallite3* relatedCrys, const Facet3* facet) 
    : m_relatedCrys((Crystallite3*)relatedCrys), facet((Facet3*)facet) {}


FrontFacet3::~FrontFacet3() {}