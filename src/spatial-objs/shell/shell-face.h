// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/edge.h"
#include "spatial-objs/vert.h"
#include "spatial-objs/shell/shell-front/shell-front-edge.h"
#include "spatial-objs/shell/shell-front/shell-front-vert.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {
namespace shell {

class Face
{
public:
    shell::Edge* edges[3];

    real_t preferredLength() const;

    const std::list<pmg::Face*>&   innerFaces() const;
    const std::list<pmg::Edge*>&   innerEdges() const;
    const std::vector<pmg::Vert*>& innerVerts() const;

    const std::list<front::Edge*>& frontEdges() const;

    shell::Vert* findVertNot( const shell::Edge* sEdge ) const;
    shell::Edge* findShellEdgeContaining( const pmg::Edge* edge ) const;

    // Needs shell edges to be already segmentized.
    void triangulate( real_t preferredLen );

    bool contains( const shell::Edge*   sEdge ) const;
    bool contains( const shell::Vert* sVert ) const;

    Face( const shell::Edge* sEdge0, const shell::Edge* sEdge1, const shell::Edge* sEdge2 );


private:
    enum class ExhaustType
    {
        DontExhaust,
        WithoutNewVert,
        WithNewVert
    };

    using FrPlEdge = front::Edge;
    using FrPlVert = front::Vert;
    using pair_dd = std::pair<real_t, real_t>;
    using pair_ff = std::pair<pmg::Face*, pmg::Face*>;
    using pair_ee = std::pair<pmg::Edge*, pmg::Edge*>;

    real_t m_preferredLength;

    std::list<pmg::Face*>  m_innerFaces;
    std::list<pmg::Edge*>   m_innerEdges;
    std::vector<pmg::Vert*> m_innerVerts;

    std::list<FrPlEdge*>   m_frontEdges;
    std::list<FrPlVert*> m_frontVerts;

    FrPlVert* findFrontVert( const pmg::Vert* vert ) const;

    FrPlEdge* addToFront( const pmg::Edge* edge );
    FrPlVert* addToFront( const pmg::Vert* vert );

    void removeFromFront( FrPlEdge* fEdge );
    void removeFromFront( FrPlVert* fVert );

    bool anyVertInsidePotentialTriangCheck( FrPlVert* fVert ) const;
    bool doesSegmentIntersectsWithFront( const Vec& p0, const Vec& p1 ) const;

    Vec computeNormalInTriang( FrPlEdge* fEdge, const Vec& oppVertPos );
    Vec computeNormalInTriang( FrPlEdge* fEdge, pmg::Edge* oneOfRemainingEdges ); // Do i need it?

    bool tryComputeNewVertPosType2( FrPlEdge* fEdge, Vec& out_pos );
    bool tryComputeNewVertPosType1( FrPlEdge* fEdge, Vec& out_pos, int smallAngleIndex );
    bool tryComputeNewVertPosType0( FrPlEdge* fEdge, Vec& out_pos );
    bool tryComputeNewVertPos(      FrPlEdge* fEdge, Vec& out_pos );

    static pair_dd computeMinMaxEdgesLengths(    const Vec& p0, const Vec& p1, const Vec& p2 );
    static pair_dd computeMinMaxEdgesSqrLengths( const Vec& p0, const Vec& p1, const Vec& p2 );
    static real_t  computeTriangSimpleQuality(     const Vec& p0, const Vec& p1, const Vec& p2 );
    static real_t  computeTriangSimpleSqrQuality(  const Vec& p0, const Vec& p1, const Vec& p2 );

    FrPlEdge* chooseEdgeForExhaustionWithNewVert( FrPlVert* fVert );
    void      exhaustWithNewVert( FrPlEdge* fEdge, const Vec& vertPos );
    void      exhaustWithoutNewVert( FrPlVert* fVert );

    bool tryExhaustWithoutNewVert( FrPlVert* fVert );
    bool tryExhaustWithNewVert(    FrPlVert* fVert );

    bool globalIntersectionCheck() const;

    FrPlVert* currentFrontVert( real_t maxCompl ) const;
    bool exhaustWithoutNewVertPriorityPredicate( FrPlVert* fEdge );
    bool exhaustWithNewVertPriorityPredicate(    FrPlVert* fEdge );
    ExhaustType computeExhaustionTypeQualityPriority(
        FrPlVert* fVert,
        FrPlEdge*& out_withNWFrontEdge, Vec*& out_withNWNewVertPos );

    void processLastFace();
    void processAngles();

    void smoothMesh( unsigned nIterations );

    pair_ff find2AdjFaces( pmg::Edge* edge ) const;
    bool flipIfNeeded( pmg::Edge* edge );
    void delaunayPostproc();

    void computeFrontNormals() const;
    void initializeFront();
};

} // namespace shell
} // namespace pmg
