// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <array>
#include <list>
#include <vector>
#include "../edge.h"
#include "../vert.h"
#include "front/edge.h"
#include "front/vert.h"
#include "../../real-type.h"
#include "../genparams.h"

#include "../../definitions.h"


namespace pmg::surface {

// TODO: separate Mesher<2>
class Face
{
public:
    // TODO: maybe use std::reference_wrapper instead of pointer
    // TODO: represent surface::Face as any number of surface::Edges instead
    std::array<surface::Edge*, 3> edges;

    real_t preferredLength() const;

    const std::list<pmg::Face*>& innerFaces() const;
    const std::list<pmg::Edge*>& innerEdges() const;
    const std::list<pmg::Vert*>& innerVerts() const;

    const std::list<front::Edge*>& frontEdges() const;

    surface::Vert* findVertNot( const surface::Edge* sEdge ) const;
    surface::Edge* findSurfaceEdgeContaining( const pmg::Edge* edge ) const;

    // Needs surface edges to be already segmentized.
    void triangulate( real_t preferredLen, genparams::Surface genParams = genparams::Surface() );
    void smoothMesh( std::size_t nIters );
    void delaunayPostP();
    void optimizeMesh( std::size_t nSmoothIters = 20, std::size_t nDelaunaySmoothIters = 3 );

    bool contains( const surface::Edge* sEdge ) const;
    bool contains( const surface::Vert* sVert ) const;

    Face( const surface::Edge* sEdge0, const surface::Edge* sEdge1, const surface::Edge* sEdge2 );


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

    real_t m_prefLen = static_cast<real_t>(0.0);

    std::list<pmg::Face*> m_innerFaces;
    std::list<pmg::Edge*> m_innerEdges;
    std::list<pmg::Vert*> m_innerVerts;

    std::list<front::Edge*> m_frontEdges;
    std::list<front::Vert*> m_frontVerts;

    front::Vert* findFrontVert( const pmg::Vert* vert ) const;

    front::Edge* addToFront( const pmg::Edge* edge );
    front::Vert* addToFront( const pmg::Vert* vert );

    void removeFromFront( front::Edge* fEdge );
    void removeFromFront( front::Vert* fVert );

    bool anyVertInsidePotentialTriangCheck( front::Vert* fVert ) const;
    bool doesSegmentIntersectsWithFront( const spt::vec3& v0, const spt::vec3& v1 ) const;
    bool doesSegmentIntersectsWithFront( const pmg::Vert* v0, const spt::vec3& v1 ) const;

    spt::vec3 computeNormalInTriang( front::Edge* fEdge, const spt::vec3& oppVertPos );
    spt::vec3 computeNormalInTriang( front::Edge* fEdge, pmg::Edge* oneOfRemainingEdges ); // Do i need it?

    bool tryComputeNewVertPosType2( front::Edge* fEdge, spt::vec3& out_pos );
    bool tryComputeNewVertPosType1( front::Edge* fEdge, spt::vec3& out_pos, std::size_t smallAngleIdx );
    bool tryComputeNewVertPosType0( front::Edge* fEdge, spt::vec3& out_pos );
    bool tryComputeNewVertPos(      front::Edge* fEdge, spt::vec3& out_pos );

    static pair_rr computeMinMaxEdgesLengths(     const spt::vec3& p0, const spt::vec3& p1, const spt::vec3& p2 );
    static pair_rr computeMinMaxEdgesSqrLengths(  const spt::vec3& p0, const spt::vec3& p1, const spt::vec3& p2 );
    static real_t  computeTriangSimpleQuality(    const spt::vec3& p0, const spt::vec3& p1, const spt::vec3& p2 );
    static real_t  computeTriangSimpleSqrQuality( const spt::vec3& p0, const spt::vec3& p1, const spt::vec3& p2 );

    front::Edge* chooseEdgeForExhaustionWithNewVert( front::Vert* fVert );
    void exhaustWithNewVert( front::Edge* fEdge, const spt::vec3& vertPos );
    void exhaustWithoutNewVert( front::Vert* fVert );

    bool tryExhaustWithoutNewVert( front::Vert* fVert );
    bool tryExhaustWithNewVert(    front::Vert* fVert );

    bool globalIntersectionCheck() const;

    front::Vert* currentFrontVert( real_t maxCompl ) const;
    bool exhaustWithoutNewVertPriorityPredicate( front::Vert* fEdge );
    bool exhaustWithNewVertPriorityPredicate(    front::Vert* fEdge );
    ExhaustType computeExhaustionTypeQualityPriority(
        front::Vert* fVert,
        front::Edge*& out_withNWFrontEdge, spt::vec3*& out_withNWNewVertPos );

    void processLastFace();
    void processAngles();

    pair_ff find2AdjFaces( pmg::Edge* edge ) const;
    bool flipIfNeeded( pmg::Edge* edge );

    void computeFrontNormals() const;
    void initializeFront();
};

} // namespace pmg::surface
