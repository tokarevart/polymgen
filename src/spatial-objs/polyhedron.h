// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <list>
#include <vector>
#include "spatial-objs/polyhedron-front/polyhedron-front-face.h"
#include "spatial-objs/polyhedron-front/polyhedron-front-edge.h"
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

    const std::list<front::Face*>& frontFaces() const;
    const std::list<front::Edge*>& frontEdges()  const;

    // Returns minimum and average tetrahedrons qualitys.
    pair_dd analyzeMeshQuality();
    void    generateMesh(real_t preferredLength);

    Polyhedron();
    Polyhedron(PolyhedralSet* polyhedr);
    ~Polyhedron();


private:
    enum class ExhaustType
    {
        DontExhaust,
        WithoutNewVert,
        WithNewVert
    };

    using pair_ff = std::pair<front::Face*, front::Face*>;

    real_t m_preferredLength;

    PolyhedralSet* m_polycr = nullptr;

    std::vector<shell::Face*> m_shellFaces;
    std::vector<shell::Edge*> m_shellEdges;
    std::vector<shell::Vert*> m_shellVerts;

    std::vector<Tetr*> m_innerTetrs;
    std::vector<Face*> m_innerFaces;
    std::vector<Edge*> m_innerEdges;
    std::vector<Vert*> m_innerVerts;

    std::list<front::Face*> m_frontFaces;
    std::list<front::Edge*> m_frontEdges;

    bool shellContains( const Vert* vert ) const;

    shell::Edge* findShellEdge( const shell::Vert* v0, const shell::Vert* v1 ) const;
    front::Face* findFrontFace( const Face* face ) const;
    std::vector<front::Edge*> findFEdge( const Vert* v0, const Vert* v1 ) const;
    std::vector<front::Edge*> findFEdge( const Edge* edge ) const;

    // Adds new front Face and corresponding Face.
    front::Face* addToFront( const Face* face, bool addInner = true );
    front::Edge* addToFront( const Edge* edge, bool addInner = true );

    void removeFromFront( front::Face* fFace );
    void removeFromFront( front::Edge* fEdge );

    bool vertInsideFrontCheck( const Vec& v ) const;
    bool segmentGlobalIntersectionCheck( const Vec& v0, const Vec& v1 ) const;
    bool segmentFrontIntersectionCheck( const Vec& v0, const Vec& v1 ) const;
    bool edgeGlobalIntersectionCheck( const Edge* edge ) const;
    bool edgeIntersectionCheck( front::Edge* fEdge ) const;
    bool faceIntersectionCheck( const Vert* v0, const Vert* v1, const Vec&    v2 ) const;
    bool faceIntersectionCheck( const Vert* v0, const Vert* v1, const Vert* v2 ) const;
    bool facesIntersectionCheck( front::Edge* fEdge ) const;
    bool insideTetrCheck( const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3, const Vec& vert ) const;
    bool anyVertInsidePotentialTetrCheck( front::Edge* fEdge ) const;
    bool parallelFacesCheck( front::Edge* fEdge ) const;
    bool doesFrontIntersectSphere( const Vec& center, real_t radius ) const;
    bool frontSplitCheck(    front::Edge* fEdge, front::Edge* oppFEdge = nullptr ) const;
    bool frontCollapseCheck( front::Edge* fEdge, front::Edge* oppFEdge = nullptr ) const;

    static pair_dd computeMinMaxEdgesLengths(    const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3 );
    static pair_dd computeMinMaxEdgesSqrLengths( const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3 );
    static real_t  computeTetrSimpleQuality(     const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3 );
    static real_t  computeTetrSimpleSqrQuality(  const Vec& p0, const Vec& p1, const Vec& p2, const Vec& p3 );

    front::Edge* currentFrontEdge( real_t maxCompl ) const;
    bool exhaustWithoutNewVertPriorityPredicate( front::Edge* fEdge );
    bool exhaustWithNewVertPriorityPredicate(    front::Edge* fEdge );
    ExhaustType computeExhaustionTypeQualityPriority(
        front::Edge* fEdge,
        front::Face*& out_withNWFFace, Vec*& out_withNWNewVertPos );

    Vec computeNormalInTetr( const front::Face* fFace, const Vec& oppVertPos ) const;
    Vec computeNormalInTetr( const front::Face* fFace, const pmg::Edge* oneOfRemainingEdges ) const;
    Vec computeNormalInTetr( const Vec& fFacePos0, const Vec& fFacePos1, const Vec& fFacePos2, const Vec& oppVertPos ) const;

    void setFEdgesInFrontSplit( const front::Edge* fEdge, front::Edge* newOppFEdges[2], front::Face* newFFaces[2], pair_ff oppFFaces ) const;
    void exhaustFrontCollapse( front::Edge* fEdge, front::Edge* oppFEdge );
    void exhaustFrontSplit( front::Edge* fEdge, front::Edge* oppFEdge );
    void exhaustWithoutNewVertOppEdgeExists( front::Edge* fEdge, front::Edge* oppFEdge );
    void exhaustWithoutNewVertOppEdgeDontExists( front::Edge* fEdge );
    void exhaustWithoutNewVert( front::Edge* fEdge, bool doesOppEdgeExists = true, front::Edge* oppFEdge = nullptr );


    bool tryComputeNewVertPosType3(
            front::Face* fFace, Vec& out_pos );
    bool tryComputeNewVertPosType2(
            front::Face* frontFace, Vec& out_pos,
            int smallAngleIndex0, int smallAngleIndex1 );
    bool tryComputeNewVertPosType1(
            front::Face* fFace, Vec& out_pos,
            int smallAngleIndex);
    bool tryComputeNewVertPosType0(
            front::Face* fFace, Vec& out_pos );
    bool tryComputeNewVertPos( front::Face* fFace, Vec& out_pos );

    real_t sqr4FaceArea( const front::Face* fFace ) const;
    front::Face* chooseFaceForExhaustionWithNewVert( front::Edge* fEdge );
    void         exhaustWithNewVert( front::Face* fFace, const Vec& vertPos );

    bool tryExhaustWithoutNewVert( front::Edge* fEdge, bool doesOppEdgeExists = true, front::Edge* oppEdge = nullptr );
    bool tryExhaustWithNewVert(    front::Edge* fEdge );

    bool globalIntersectionCheck();

    bool isFrontExhausted();
    void processAngles();

    void debug();

    void smoothMesh(           unsigned nIterations );
    void smoothNotFinisedMesh( unsigned nIterations );
    void smoothFront(          unsigned nIterations );
    void smoothAroundFrontVert( Vert* frontVert );

    void computeFrontNormals();
    void initializeFFaceFEdges( front::Face* fFace ) const;
    void initializeFront();
};

} // namespace pmg
