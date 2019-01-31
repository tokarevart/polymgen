#pragma once
#include "Crystallite3.h"
#include "Edge3.h"
#include "Vertex3.h"
#include "Helpers/SpatialAlgs/Vec3.h"

class Crystallite3;
class Edge3;
class Vertex3;
namespace tva { struct Vec3; }


class Facet3
{
public:
    Edge3* edges[3];

    tva::Vec3 computeCenter();
    double computeQuality();

    static Edge3* intersectAlongAnEdge(const Facet3* facet0, const Facet3* facet1);

    bool intersectsBy(const tva::Point3& origin, const tva::Vec3& dir);
    Vertex3* findVertNotIncludedInEdge(const Edge3* edge) const;
    Edge3*   findEdgeNotContainingVert(const Vertex3* vert) const;
    Edge3*   findEdge(const Vertex3* vert0, const Vertex3* vert1);
    Edge3* findShortestEdge();
    Edge3* findLongestEdge();

    bool contains(const Edge3*   edge) const;
    bool contains(const Vertex3* vert) const;

    Facet3();
    Facet3(const Edge3* edge0, const Edge3* edge1, const Edge3* edge2);
    Facet3(const Vertex3* vert0, const Vertex3* vert1, const Vertex3* vert2);
    ~Facet3();
};


class FrontFacet3
{
public:
    Facet3* facet;

    tva::Vec3 getNormal();
    void      setNormal(const tva::Vec3& vec);
    tva::Vec3 computeNormal();

    tva::Vec3 computeCenter();
    double    computeQuality();

    FrontFacet3(const Crystallite3* relatedCrys, const Facet3* facet);
    ~FrontFacet3();


private:
    Crystallite3* m_relatedCrys;
    tva::Vec3 m_normal;
};