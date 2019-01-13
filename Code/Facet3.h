#pragma once
#include <list>
#include <memory>
#include "Definitions.h"
#include "Inclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;

class Facet3 : public unique_ptr_helper<Facet3>
{
public:
    unique_ptr<Edge3>* edges[3];

    const Vec3 computeCenter();
    const double computeQuality();

    static unique_ptr<Edge3>* intersectAlongAnEdge(const Facet3& facet0, const Facet3& facet1);

    const bool intersectsBy(const Point3& origin, const Vec3& dir);
    unique_ptr<Vertex3>* findVertexNotIncludedInEdge(const Edge3& edge) const;
    unique_ptr<Edge3>*   findEdgeNotContainingVertex(const Vertex3& vert) const;
    unique_ptr<Edge3>*   findEdge(const Vertex3& vert0, const Vertex3& vert1);
    unique_ptr<Edge3>* findShortestEdge();
    unique_ptr<Edge3>* findLongestEdge();

    const bool contains(const Edge3&   edge) const;
    const bool contains(const Vertex3& vert) const;

    Facet3();
    Facet3(Edge3& edge0,   Edge3& edge1,   Edge3& edge2);
    Facet3(Vertex3& vert0, Vertex3& vert1, Vertex3& vert2);
    ~Facet3();
};

class FrontFacet3 : public unique_ptr_helper<FrontFacet3>
{
    Vec3 _normal;

public:
    Facet3* facet;

    const Vec3 getNormal();
    const void setNormal(const Vec3& vec);
    const Vec3 computeNormal(const vector<unique_ptr<FrontFacet3>*>& frontFacets);

    const Vec3 computeCenter();
    const double computeQuality();

    FrontFacet3(Facet3* facet);
    ~FrontFacet3();
};