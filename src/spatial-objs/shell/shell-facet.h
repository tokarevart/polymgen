#pragma once
#include <list>
#include <vector>
#include "spatial-objs/edge.h"
#include "spatial-objs/vertex.h"
#include "spatial-objs/front/plane/front-plane-edge.h"
#include "spatial-objs/front/plane/front-plane-vertex.h"

#include "definitions.h"


namespace pmg {
namespace shell {

class Facet
{
public:
    shell::Edge* edges[3];

    double preferredLength() const;

    const std::list<pmg::Facet*>&  innerFacets() const;
    const std::list<pmg::Edge*>&   innerEdges()  const;
    const std::vector<pmg::Vertex*>& innerVerts()  const;

    const std::list<front::plane::Edge*>& frontEdges() const;

    shell::Vertex* findVertNot( const shell::Edge* sEdge ) const;
    shell::Edge*   findShellEdgeContaining(   const   pmg::Edge*  edge ) const;

    // Needs shell edges to be already segmentized.
    void triangulate( double preferredLen );

    bool contains( const shell::Edge*   sEdge ) const;
    bool contains( const shell::Vertex* sVert ) const;

    Facet( const shell::Edge* sEdge0, const shell::Edge* sEdge1, const shell::Edge* sEdge2 );


private:
    enum class ExhaustType
    {
        DontExhaust,
        WithoutNewVert,
        WithNewVert
    };

    using FrPlEdge   = front::plane::Edge;
    using FrPlVertex = front::plane::Vertex;
    using pair_dd = std::pair<double, double>;
    using pair_ff = std::pair<pmg::Facet*, pmg::Facet*>;
    using pair_ee = std::pair<pmg::Edge*, pmg::Edge*>;
    using Vec   = tva::Vec;
    using Point = tva::Point;

    double m_preferredLength;

    std::list<pmg::Facet*>  m_innerFacets;
    std::list<pmg::Edge*>   m_innerEdges;
    std::vector<pmg::Vertex*> m_innerVerts;

    std::list<FrPlEdge*>   m_frontEdges;
    std::list<FrPlVertex*> m_frontVerts;

    FrPlVertex* findFrontVert( const pmg::Vertex* vert ) const;

    FrPlEdge*   addToFront( const pmg::Edge* edge );
    FrPlVertex* addToFront( const pmg::Vertex* vert );

    void removeFromFront( FrPlEdge* fEdge );
    void removeFromFront( FrPlVertex* fVert );

    bool anyVertInsidePotentialTriangCheck( FrPlVertex* fVert ) const;
    bool doesSegmentIntersectsWithFront( const Point& p0, const Point& p1 ) const;

    Vec computeNormalInTriang( FrPlEdge* fEdge, const Point& oppVertPos );
    Vec computeNormalInTriang( FrPlEdge* fEdge, pmg::Edge* oneOfRemainingEdges ); // Do i need it?

    bool tryComputeNewVertPosType2( FrPlEdge* fEdge, Point& out_pos );
    bool tryComputeNewVertPosType1( FrPlEdge* fEdge, Point& out_pos, int smallAngleIndex );
    bool tryComputeNewVertPosType0( FrPlEdge* fEdge, Point& out_pos );
    bool tryComputeNewVertPos(      FrPlEdge* fEdge, Point& out_pos );

    static pair_dd computeMinMaxEdgesLengths(    const Vec& p0, const Vec& p1, const Vec& p2 );
    static pair_dd computeMinMaxEdgesSqrLengths( const Vec& p0, const Vec& p1, const Vec& p2 );
    static double  computeTriangSimpleQuality(     const Vec& p0, const Vec& p1, const Vec& p2 );
    static double  computeTriangSimpleSqrQuality(  const Vec& p0, const Vec& p1, const Vec& p2 );

    FrPlEdge* chooseEdgeForExhaustionWithNewVert( FrPlVertex* fVert );
    void      exhaustWithNewVert( FrPlEdge* fEdge, const Point& vertPos );
    void      exhaustWithoutNewVert( FrPlVertex* fVert );

    bool tryExhaustWithoutNewVert( FrPlVertex* fVert );
    bool tryExhaustWithNewVert(    FrPlVertex* fVert );

    bool globalIntersectionCheck() const;

    FrPlVertex* currentFrontVert( double maxCompl ) const;
    bool exhaustWithoutNewVertPriorityPredicate( FrPlVertex* fEdge );
    bool exhaustWithNewVertPriorityPredicate(    FrPlVertex* fEdge );
    ExhaustType computeExhaustionTypeQualityPriority(
        FrPlVertex* fVert,
        FrPlEdge*& out_withNWFrontEdge, Point*& out_withNWNewVertPos );

    void processLastFacet();
    void processAngles();

    void smoothMesh( unsigned nIterations );

    pair_ff find2AdjFacets( pmg::Edge* edge ) const;
    bool flipIfNeeded( pmg::Edge* edge );
    void delaunayPostproc();

    void computeFrontNormals() const;
    void initializeFront();
};

} // namespace shell
} // namespace pmg
