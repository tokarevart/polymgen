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
#include "../pmg-settings.h"

#include "definitions.h"


namespace pmg {
namespace shell {

class Face
{
public:
    shell::Edge* edges[3];

    real_t preferredLength() const;

    const std::list  <pmg::Face*>& innerFaces() const;
    const std::list  <pmg::Edge*>& innerEdges() const;
    const std::vector<pmg::Vert*>& innerVerts() const;

    const std::list<front::Edge*>& frontEdges() const;

    shell::Vert* findVertNot( const shell::Edge* sEdge ) const;
    shell::Edge* findShellEdgeContaining( const pmg::Edge* edge ) const;

    // Needs shell edges to be already segmentized.
    void triangulate( real_t preferredLen, gensettings::Shell genSettings = gensettings::Shell() );
    void smoothMesh( size_t nIters );
    void delaunayPostP();
    void optimizeMesh( size_t nSmoothIters = 20, size_t nDelaunaySmoothIters = 3 );

    bool contains( const shell::Edge* sEdge ) const;
    bool contains( const shell::Vert* sVert ) const;

    Face( const shell::Edge* sEdge0, const shell::Edge* sEdge1, const shell::Edge* sEdge2 );


private:
    enum class ExhaustType
    {
        DontExhaust,
        WithoutNewVert,
        WithNewVert
    };

    using pair_rr = std::pair<real_t, real_t>;
    using pair_ff = std::pair<pmg::Face*, pmg::Face*>;
    using pair_ee = std::pair<pmg::Edge*, pmg::Edge*>;

    real_t m_prefLen;

    std::list<pmg::Face*>   m_innerFaces;
    std::list<pmg::Edge*>   m_innerEdges;
    std::vector<pmg::Vert*> m_innerVerts;

    std::list<front::Edge*> m_frontEdges;
    std::list<front::Vert*> m_frontVerts;

    front::Vert* findFrontVert( const pmg::Vert* vert ) const;

    front::Edge* addToFront( const pmg::Edge* edge );
    front::Vert* addToFront( const pmg::Vert* vert );

    void removeFromFront( front::Edge* fEdge );
    void removeFromFront( front::Vert* fVert );

    bool anyVertInsidePotentialTriangCheck( front::Vert* fVert ) const;
    bool doesSegmentIntersectsWithFront( const vec3& p0, const vec3& p1 ) const;

    vec3 computeNormalInTriang( front::Edge* fEdge, const vec3& oppVertPos );
    vec3 computeNormalInTriang( front::Edge* fEdge, pmg::Edge* oneOfRemainingEdges ); // Do i need it?

    bool tryComputeNewVertPosType2( front::Edge* fEdge, vec3& out_pos );
    bool tryComputeNewVertPosType1( front::Edge* fEdge, vec3& out_pos, int smallAngleIndex );
    bool tryComputeNewVertPosType0( front::Edge* fEdge, vec3& out_pos );
    bool tryComputeNewVertPos(      front::Edge* fEdge, vec3& out_pos );

    static pair_rr computeMinMaxEdgesLengths(     const vec3& p0, const vec3& p1, const vec3& p2 );
    static pair_rr computeMinMaxEdgesSqrLengths(  const vec3& p0, const vec3& p1, const vec3& p2 );
    static real_t  computeTriangSimpleQuality(    const vec3& p0, const vec3& p1, const vec3& p2 );
    static real_t  computeTriangSimpleSqrQuality( const vec3& p0, const vec3& p1, const vec3& p2 );

    front::Edge* chooseEdgeForExhaustionWithNewVert( front::Vert* fVert );
    void exhaustWithNewVert( front::Edge* fEdge, const vec3& vertPos );
    void exhaustWithoutNewVert( front::Vert* fVert );

    bool tryExhaustWithoutNewVert( front::Vert* fVert );
    bool tryExhaustWithNewVert(    front::Vert* fVert );

    bool globalIntersectionCheck() const;

    front::Vert* currentFrontVert( real_t maxCompl ) const;
    bool exhaustWithoutNewVertPriorityPredicate( front::Vert* fEdge );
    bool exhaustWithNewVertPriorityPredicate(    front::Vert* fEdge );
    ExhaustType computeExhaustionTypeQualityPriority(
        front::Vert* fVert,
        front::Edge*& out_withNWFrontEdge, vec3*& out_withNWNewVertPos );

    void processLastFace();
    void processAngles();

    pair_ff find2AdjFaces( pmg::Edge* edge ) const;
    bool flipIfNeeded( pmg::Edge* edge );

    void computeFrontNormals() const;
    void initializeFront();
};

} // namespace shell
} // namespace pmg
