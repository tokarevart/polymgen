// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/front/surface/front-surface-facet.h"
#include "spatial-objs/front/surface/front-surface-edge.h"
#include "spatial-objs/shell/shell-facet.h"
#include "spatial-objs/shell/shell-edge.h"
#include "spatial-objs/shell/shell-vertex.h"
#include "spatial-objs/tetr.h"
#include "spatial-objs/facet.h"
#include "spatial-objs/edge.h"
#include "spatial-objs/vertex.h"
#include "helpers/spatial-algs/vec.h"

#include "spatial-objs/polycrystal.h"

#include "definitions.h"


namespace pmg {

class Crystallite
{
    using FrSuFacet = front::surface::Facet;
    using FrSuEdge  = front::surface::Edge;
    using pair_dd = std::pair<double, double>;

public:
    double preferredLength() const;

    void addToShell(const shell::Facet*  shellFacet);
    void addToShell(const shell::Edge*   shellEdge);
    void addToShell(const shell::Vertex* shellVert);

    bool shellContains(const shell::Facet*  shellFacet) const;
    bool shellContains(const shell::Edge*   shellEdge)  const;
    bool shellContains(const shell::Vertex* shellVert)  const;

    const std::vector<Tetr*>&   innerTetrs()  const;
    const std::vector<Facet*>&  innerFacets() const;
    const std::vector<Vertex*>& innerVerts()  const;

    const std::list<FrSuFacet*>& frontFacets() const;
    const std::list<FrSuEdge*>&  frontEdges()  const;

    // Returns minimum and average tetrahedrons qualitys.
    pair_dd analyzeMeshQuality();
    void    generateMesh(double preferredLength);

    Crystallite();
    Crystallite(Polycrystal* polycr);
    ~Crystallite();


private:
    enum class ExhaustType
    {
        DontExhaust,
        WithoutNewVert,
        WithNewVert
    };

    using pair_ff = std::pair<FrSuFacet*, FrSuFacet*>;
    using Vec   = tva::Vec;
    using Point = tva::Point;

    double m_preferredLength;

    Polycrystal* m_polycr = nullptr;

    std::vector<shell::Facet*>  m_shellFacets;
    std::vector<shell::Edge*>   m_shellEdges;
    std::vector<shell::Vertex*> m_shellVerts;

    std::vector<Tetr*>    m_innerTetrs;
    std::vector<Facet*>   m_innerFacets;
    std::vector<Edge*>    m_innerEdges;
    std::vector<Vertex*>  m_innerVerts;

    std::list<FrSuFacet*> m_frontFacets;
    std::list<FrSuEdge*>  m_frontEdges;

    bool shellContains( const Vertex* vert ) const;

    shell::Edge* findShellEdge( const shell::Vertex* v0, const shell::Vertex* v1 ) const;
    FrSuFacet* findFrontFacet( const Facet* facet ) const;
    std::vector<FrSuEdge*>  findFEdge( const Vertex* v0, const Vertex* v1 ) const;
    std::vector<FrSuEdge*>  findFEdge( const Edge* edge ) const;

    // Adds new front facet and corresponding facet.
    FrSuFacet* addToFront( const Facet* facet, bool addInner = true );
    FrSuEdge*  addToFront( const Edge*  edge,  bool addInner = true );

    void removeFromFront( FrSuFacet* fFacet );
    void removeFromFront( FrSuEdge* fEdge );

    bool vertInsideFrontCheck( const Point& v ) const;
    bool segmentGlobalIntersectionCheck( const Point& v0, const Point& v1 ) const;
    bool segmentFrontIntersectionCheck( const Point& v0, const Point& v1 ) const;
    bool edgeGlobalIntersectionCheck( const Edge* edge ) const;
    bool edgeIntersectionCheck( FrSuEdge* fEdge ) const;
    bool facetIntersectionCheck( const Vertex* v0, const Vertex* v1, const Vec&    v2 ) const;
    bool facetIntersectionCheck( const Vertex* v0, const Vertex* v1, const Vertex* v2 ) const;
    bool facetsIntersectionCheck( FrSuEdge* fEdge ) const;
    bool insideTetrCheck( const Point& p0, const Point& p1, const Point& p2, const Point& p3, const Point& vert ) const;
    bool anyVertInsidePotentialTetrCheck( FrSuEdge* fEdge ) const;
    bool parallelFacetsCheck( FrSuEdge* fEdge ) const;
    bool doesFrontIntersectSphere( const Point& center, double radius ) const;
    bool frontSplitCheck(    FrSuEdge* fEdge, FrSuEdge* oppFEdge = nullptr ) const;
    bool frontCollapseCheck( FrSuEdge* fEdge, FrSuEdge* oppFEdge = nullptr ) const;

    static pair_dd computeMinMaxEdgesLengths(    const Point& p0, const Point& p1, const Point& p2, const Point& p3 );
    static pair_dd computeMinMaxEdgesSqrLengths( const Point& p0, const Point& p1, const Point& p2, const Point& p3 );
    static double  computeTetrSimpleQuality(     const Point& p0, const Point& p1, const Point& p2, const Point& p3 );
    static double  computeTetrSimpleSqrQuality(  const Point& p0, const Point& p1, const Point& p2, const Point& p3 );

    FrSuEdge* currentFrontEdge( double maxCompl ) const;
    bool exhaustWithoutNewVertPriorityPredicate( FrSuEdge* fEdge );
    bool exhaustWithNewVertPriorityPredicate(    FrSuEdge* fEdge );
    ExhaustType computeExhaustionTypeQualityPriority(
        FrSuEdge* fEdge,
        FrSuFacet*& out_withNWFFacet, Vec*& out_withNWNewVertPos );

    Vec computeNormalInTetr( const FrSuFacet* fFacet, const Point& oppVertPos ) const;
    Vec computeNormalInTetr( const FrSuFacet* fFacet, const pmg::Edge* oneOfRemainingEdges ) const;
    Vec computeNormalInTetr( const Point& fFacetPos0, const Point& fFacetPos1, const Point& fFacetPos2, const Point& oppVertPos ) const;

    void setFEdgesInFrontSplit( const FrSuEdge* fEdge, FrSuEdge* newOppFEdges[2], FrSuFacet* newFFacets[2], pair_ff oppFFacets ) const;
    void exhaustFrontCollapse( FrSuEdge* fEdge, FrSuEdge* oppFEdge );
    void exhaustFrontSplit( FrSuEdge* fEdge, FrSuEdge* oppFEdge );
    void exhaustWithoutNewVertOppEdgeExists( FrSuEdge* fEdge, FrSuEdge* oppFEdge );
    void exhaustWithoutNewVertOppEdgeDontExists( FrSuEdge* fEdge );
    void exhaustWithoutNewVert( FrSuEdge* fEdge, bool doesOppEdgeExists = true, FrSuEdge* oppFEdge = nullptr );


    bool tryComputeNewVertPosType3(
            FrSuFacet* fFacet, Vec& out_pos );
    bool tryComputeNewVertPosType2(
            FrSuFacet* frontFacet, Vec& out_pos,
            int smallAngleIndex0, int smallAngleIndex1 );
    bool tryComputeNewVertPosType1(
            FrSuFacet* fFacet, Vec& out_pos,
            int smallAngleIndex);
    bool tryComputeNewVertPosType0(
            FrSuFacet* fFacet, Vec& out_pos );
    bool tryComputeNewVertPos( FrSuFacet* fFacet, Vec& out_pos );

    double sqr4FacetArea( const FrSuFacet* fFacet ) const;
    // Сurrently i don't know which method of choice is better
    FrSuFacet* chooseFacetForExhaustionWithNewVert( FrSuEdge* fEdge );
    void       exhaustWithNewVert( FrSuFacet* fFacet, const Vec& vertPos );

    bool tryExhaustWithoutNewVert( FrSuEdge* fEdge, bool doesOppEdgeExists = true, FrSuEdge* oppEdge = nullptr );
    bool tryExhaustWithNewVert(    FrSuEdge* fEdge ); // Note: check optimization

    bool globalIntersectionCheck();

    bool isFrontExhausted();
    void processAngles();

    void debug();

    void smoothMesh(           unsigned nIterations );
    void smoothNotFinisedMesh( unsigned nIterations );
    void smoothFront(          unsigned nIterations );
    void smoothAroundFrontVert( Vertex* frontVert );

    void computeFrontNormals();
    void initializeFFacetFEdges( FrSuFacet* fFacet ) const;
    void initializeFront();
};

} // namespace pmg
