#pragma once
#include <list>
#include <vector>
#include <memory>
#include "Definitions.h"
#include "Inclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;
using std::list;
using std::vector;

class Edge3 : public unique_ptr_helper<Edge3>
{
public:
    unique_ptr<Vertex3>* vertexes[2];

    const double magnitude() const;
    const double sqrMagnitude() const;

    void flip(
        vector<unique_ptr<Edge3>*>& edges, 
        vector<unique_ptr<Facet3>*>& facets);
    const bool flipIfNeeded(
        vector<unique_ptr<Edge3>*>& edges, 
        vector<unique_ptr<Facet3>*>& facets);
    // Adj means adjacent.
    void findAdjFacets(
        const vector<unique_ptr<Facet3>*>& facets,
        vector<unique_ptr<Facet3>*>& adjFacets);
    void find2AdjFacets(
        const vector<unique_ptr<Facet3>*>& facets, 
        unique_ptr<Facet3>*& facet0, 
        unique_ptr<Facet3>*& facet1);
    // Doesn't change vertexes relations.
    void make2Instead(
        vector<unique_ptr<Facet3>*>& facets, 
        vector<unique_ptr<Edge3>*>& edges, 
        vector<unique_ptr<Vertex3>*>& vertexes);

    const bool contains(const Vertex3& vertex) const;

    const bool belongsToShell();

    const bool needToFlip(const vector<unique_ptr<Facet3>*>& facets);
    
    Edge3();
    Edge3(Vertex3& vertex0, Vertex3& vertex1);
    ~Edge3();
};

class FrontEdge3 : public unique_ptr_helper<FrontEdge3>
{
    double _exCos;
    // For the future.
    //unique_ptr<FrontFacet3>* _adjFacets[2];

public:
    Edge3* edge;
    bool needProcessing = true;

    double getAngleExCos    (const Crystallite3* crys);
    double computeAngleExCos(const Crystallite3* crys);
    double computeAngleExCos_OLD(const Crystallite3* crys);
    double computeAngleCos  (bool &out_isConcave, const Crystallite3* crys);
    double computeAngle     (const Crystallite3* crys);

    // Adj means adjacent.
    const bool adjFacetsArrayContains(const unique_ptr<FrontFacet3>* frontFacet);
    unique_ptr<FrontFacet3>** getAdjFacets();
    const bool addAdjFacetInArray(const unique_ptr<FrontFacet3>* frontFacet);
    const bool removeAdjFacetFromArray(const unique_ptr<FrontFacet3>* frontFacet);

    //void findAdjacentFrontFacets(
    //    const vector<unique_ptr<FrontFacet3>*>& frontFacets, 
    //    unique_ptr<FrontFacet3>*& out_frontFacet0,
    //    unique_ptr<FrontFacet3>*& out_frontFacet1);
    //void findOppositeVertexes(
    //    const vector<unique_ptr<FrontFacet3>*>& frontFacets, 
    //    const vector<unique_ptr<FrontEdge3>*>&  frontEdges, 
    //    unique_ptr<Vertex3>*& out_vert0, unique_ptr<Vertex3>*& out_vert1);
    //unique_ptr<FrontEdge3>* findOppositeFrontEdge(
    //    const vector<unique_ptr<FrontFacet3>*>& frontFacets,
    //    const vector<unique_ptr<FrontEdge3>*>& frontEdges);

    FrontEdge3(Edge3* edge);
    ~FrontEdge3();
};