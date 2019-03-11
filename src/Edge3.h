#pragma once
#include <list>
#include <vector>
#include "Crystallite3.h"
#include "Facet3.h"
#include "Vertex3.h"

class Crystallite3;
class Facet3;
class FrontFacet3;
class Vertex3;


class Edge3
{
public:
    Vertex3* verts[2];

    double magnitude()    const;
    double sqrMagnitude() const;

    void flip(
        std::list<Edge3*>& edgesList,
        std::list<Facet3*>& facetsList);
    bool flipIfNeeded(
        std::list<Edge3*>& edgesList,
        std::list<Facet3*>& facetsList);
    // Adj means adjacent.
    void findAdjFacets(
        const std::list<Facet3*>& facetsList,
        std::list<Facet3*>& adjFacets);
    void find2AdjFacets(
        const std::list<Facet3*>& facetsList,
        Facet3*& out_facet0, 
        Facet3*& out_facet1);
    // Doesn't change vertexes relations.
    void make2Instead(
        std::list<Facet3*>& facetsList,
        std::list<Edge3*>& edgesList,
        std::list<Vertex3*>& vertsList);

    bool contains(const Vertex3* vert) const;

    bool belongsToShell();

    bool needToFlip(const std::list<Facet3*>& facetsList);
    
    Edge3();
    Edge3(const Vertex3* vert0, const Vertex3* vert1);
    ~Edge3();
};


class FrontEdge3
{
public:
    Edge3* edge;

    void refreshAngleData();

    double getComplexity();
    double getAngleExCos();
    double computeComplexity();
    double computeAngleExCos();
    double computeAngle();

    // Opp means opposite.
    std::pair<Vertex3*, Vertex3*> findOppVerts();
    FrontEdge3* findOppEdge();
    
    // Adj means adjacent.
    std::pair<FrontFacet3*, FrontFacet3*> getAdjFacets();

    bool addAdjFacetInPair(const FrontFacet3* fFacet);
    bool removeAdjFacetFromPair(const FrontFacet3* fFacet);
    bool adjFacetsPairContains(const FrontFacet3* fFacet) const;
    
    FrontEdge3(const Crystallite3* relatedCrys, const Edge3* edge);
    ~FrontEdge3();


private:
    Crystallite3* m_relatedCrys;
    std::pair<FrontFacet3*, FrontFacet3*> m_adjFacets{ nullptr, nullptr };

    bool m_needAngleExCosProcessing = true;
    bool m_needComplexityProcessing = true;
    double m_exCos;
    double m_complexity;

    bool isAdjFacetsFull();
    std::pair<FrontFacet3*, FrontFacet3*> fillAdjFacets();
};