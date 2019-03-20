// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/front/surface/front-surface-face.h"
#include "spatial-objs/front/surface/front-surface-edge.h"
#include "spatial-objs/shell/shell-face.h"
#include "spatial-objs/shell/shell-edge.h"
#include "spatial-objs/shell/shell-vertex.h"
#include "spatial-objs/tetr.h"
#include "spatial-objs/face.h"
#include "spatial-objs/edge.h"
#include "spatial-objs/vert.h"
#include "helpers/spatial-algs/vec.h"
#include "real-type.h"

#include "spatial-objs/polyhedral-set.h"

#include "definitions.h"


namespace pmg {

class Polyhedron
{
    using FrSuFace = front::surface::Face;
    using FrSuEdge  = front::surface::Edge;
    using pair_dd = std::pair<real_t, real_t>;

public:
    real_t preferredLength() const;

    void addToShell(const shell::Face* shellFace);
    void addToShell(const shell::Edge* shellEdge);
    void addToShell(const shell::Vert* shellVert);

    bool shellContains(const shell::Face* shellFace) const;
    bool shellContains(const shell::Edge* shellEdge)  const;
    bool shellContains(const shell::Vert* shellVert)  const;

    const std::vector<Tetr*>& innerTetrs()  const;
    const std::vector<Face*>& innerFaces() const;
    const std::vector<Vert*>& innerVerts()  const;

    const std::list<FrSuFace*>& frontFaces() const;
    const std::list<FrSuEdge*>& frontEdges()  const;

    // Returns minimum and average tetrahedrons qualitys.
    pair_dd analyzeMeshQuality();
    void    generateMesh(real_t preferredLength);

    Polyhedron();
    Polyhedron(PolyhedralSet* polycr);
    ~Polyhedron();


private:
    enum class ExhaustType
    {
        DontExhaust,
        WithoutNewVert,
        WithNewVert
    };

    using pair_ff = std::pair<FrSuFace*, FrSuFace*>;

    real_t m_preferredLength;

    PolyhedralSet* m_polycr = nullptr;

    std::vector<shell::Face*> m_shellFaces;
    std::vector<shell::Edge*> m_shellEdges;
    std::vector<shell::Vert*> m_shellVerts;

    std::vector<Tetr*> m_innerTetrs;
    std::vector<Face*> m_innerFaces;
    std::vector<Edge*> m_innerEdges;
    std::vector<Vert*> m_innerVerts;

    std::list<FrSuFace*> m_frontFaces;
    std::list<FrSuEdge*> m_frontEdges;

    bool shellContains( const Vert* vert ) const;

    shell::Edge* findShellEdge( const shell::Vert* v0, const shell::Vert* v1 ) const;
    FrSuFace* findFrontFace( const Face* face ) const;
    std::vector<FrSuEdge*> findFEdge( const Vert* v0, const Vert* v1 ) const;
    std::vector<FrSuEdge*> findFEdge( const Edge* edge ) const;

    // Adds new front Face and corresponding Face.
    FrSuFace* addToFront( const Face* face, bool addInner = true );
    FrSuEdge* addToFront( const Edge* edge, bool addInner = true );

    void removeFromFront( FrSuFace* fFace );
    void removeFromFront( FrSuEdge* fEdge );

    bool vertInsideFrontCheck( const Vec& v ) const;
    bool segmentGlobalIntersectionCheck( const Vec& v0, const Vec& v1 ) const;
    bool segmentFrontIntersectionCheck( const Vec& v0, const Vec& v1 ) const;
    bool edgeGlobalIntersectionCheck( const Edge* edge ) const;
    bool edgeIntersectionCheck( FrSuEdge* fEdge ) const;
    bool faceIntersectionCheck( const Vert* v0, const Vert* v1, const Vec&    v2 ) const;
    bool faceIntersectionCheck( const Vert* v0, const Vert* v1, const Vert* v2 ) const;
    bool facesIntersectionCheck( FrSuEdge* fEdge ) const;
    bool insideTetrCheck( const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3, const Vec& vert ) const;
    bool anyVertInsidePotentialTetrCheck( FrSuEdge* fEdge ) const;
    bool parallelFacesCheck( FrSuEdge* fEdge ) const;
    bool doesFrontIntersectSphere( const Vec& center, real_t radius ) const;
    bool frontSplitCheck(    FrSuEdge* fEdge, FrSuEdge* oppFEdge = nullptr ) const;
    bool frontCollapseCheck( FrSuEdge* fEdge, FrSuEdge* oppFEdge = nullptr ) const;

    static pair_dd computeMinMaxEdgesLengths(    const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3 );
    static pair_dd computeMinMaxEdgesSqrLengths( const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3 );
    static real_t  computeTetrSimpleQuality(     const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3 );
    static real_t  computeTetrSimpleSqrQuality(  const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3 );

    FrSuEdge* currentFrontEdge( real_t maxCompl ) const;
    bool exhaustWithoutNewVertPriorityPredicate( FrSuEdge* fEdge );
    bool exhaustWithNewVertPriorityPredicate(    FrSuEdge* fEdge );
    ExhaustType computeExhaustionTypeQualityPriority(
        FrSuEdge* fEdge,
        FrSuFace*& out_withNWFFace, Vec*& out_withNWNewVertPos );

    Vec computeNormalInTetr( const FrSuFace* fFace, const Vec& oppVertPos ) const;
    Vec computeNormalInTetr( const FrSuFace* fFace, const pmg::Edge* oneOfRemainingEdges ) const;
    Vec computeNormalInTetr( const Vec& fFacePos0, const Vec& fFacePos1, const Vec& fFacePos2, const Vec& oppVertPos ) const;

    void setFEdgesInFrontSplit( const FrSuEdge* fEdge, FrSuEdge* newOppFEdges[2], FrSuFace* newFFaces[2], pair_ff oppFFaces ) const;
    void exhaustFrontCollapse( FrSuEdge* fEdge, FrSuEdge* oppFEdge );
    void exhaustFrontSplit( FrSuEdge* fEdge, FrSuEdge* oppFEdge );
    void exhaustWithoutNewVertOppEdgeExists( FrSuEdge* fEdge, FrSuEdge* oppFEdge );
    void exhaustWithoutNewVertOppEdgeDontExists( FrSuEdge* fEdge );
    void exhaustWithoutNewVert( FrSuEdge* fEdge, bool doesOppEdgeExists = true, FrSuEdge* oppFEdge = nullptr );


    bool tryComputeNewVertPosType3(
            FrSuFace* fFace, Vec& out_pos );
    bool tryComputeNewVertPosType2(
            FrSuFace* frontFace, Vec& out_pos,
            int smallAngleIndex0, int smallAngleIndex1 );
    bool tryComputeNewVertPosType1(
            FrSuFace* fFace, Vec& out_pos,
            int smallAngleIndex);
    bool tryComputeNewVertPosType0(
            FrSuFace* fFace, Vec& out_pos );
    bool tryComputeNewVertPos( FrSuFace* fFace, Vec& out_pos );

    real_t sqr4FaceArea( const FrSuFace* fFace ) const;
    // Сurrently i don't know which method of choice is better
    FrSuFace* chooseFaceForExhaustionWithNewVert( FrSuEdge* fEdge );
    void       exhaustWithNewVert( FrSuFace* fFace, const Vec& vertPos );

    bool tryExhaustWithoutNewVert( FrSuEdge* fEdge, bool doesOppEdgeExists = true, FrSuEdge* oppEdge = nullptr );
    bool tryExhaustWithNewVert(    FrSuEdge* fEdge ); // Note: check optimization

    bool globalIntersectionCheck();

    bool isFrontExhausted();
    void processAngles();

    void debug();

    void smoothMesh(           unsigned nIterations );
    void smoothNotFinisedMesh( unsigned nIterations );
    void smoothFront(          unsigned nIterations );
    void smoothAroundFrontVert( Vert* frontVert );

    void computeFrontNormals();
    void initializeFFaceFEdges( FrSuFace* fFace ) const;
    void initializeFront();
};

} // namespace pmg
