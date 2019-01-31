#pragma once
#include <list>
#include <vector>
#include "Simplex3.h"
#include "Facet3.h"
#include "ShellFacet3.h"
#include "Edge3.h"
#include "ShellEdge3.h"
#include "Vertex3.h"
#include "ShellVertex3.h"
#include "Helpers/SpatialAlgs/Vec3.h"

class Simplex3;
class FrontFacet3;
class ShellFacet3;
class Facet3;
class FrontEdge3;
class ShellEdge3;
class Edge3;
class Vertex3;
class ShellVertex3;
namespace tva { struct Vec3; }

using tva::Vec3;


class Crystallite3
{
public:
    double getPreferredLength();

    void addShellFacet(const ShellFacet3* shellFacet);
    void addShellEdge(const ShellEdge3* shellEdge);

    bool shellFacetsContains(const ShellFacet3* shellFacet);
    bool shellEdgesContains(const ShellEdge3* shellEdge);

    const std::vector<Simplex3*>& getInnerSimplexes3();
    const std::vector<Facet3*>& getInnerFacets();
    const std::vector<Vertex3*>& getInnerVertexes();

    const std::list<FrontFacet3*>& getFrontFacets();
    const std::list<FrontEdge3*>& getFrontEdges();

    std::pair<double, double> analyzeMeshQuality();
    void generateMesh(double preferredLength);
    void computeFrontNormals();
    void setStartFront(const std::list<Edge3*>& edges, const std::list<Facet3*>& facets);

    Crystallite3();
    ~Crystallite3();


private:
    enum ExhaustType
    {
        DONT_EXHAUST,
        WITHOUT_NEW_VERTEX,
        WITH_NEW_VERTEX
    };

    double m_preferredLength;

    std::vector<ShellFacet3*> m_shellFacets;
    std::vector<ShellEdge3*>  m_shellEdges;

    std::vector<Simplex3*> m_innerSimps;
    std::vector<Facet3*>   m_innerFacets;
    std::vector<Edge3*>    m_innerEdges;
    std::vector<Vertex3*>  m_innerVerts;

    std::list<FrontFacet3*> m_frontFacets;
    std::list<FrontEdge3*>  m_frontEdges;

    ShellEdge3* findShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const;

    FrontFacet3* findFrontFacet(const Facet3* facet);
    FrontEdge3* findFrontEdge(const Vertex3* v0, const Vertex3* v1);
    FrontEdge3* findFrontEdge(const Edge3* edge);

    // Adds new front facet and corresponding facet.
    FrontFacet3* addFrontFacet3(const Edge3* edge0, const Edge3* edge1, const Edge3* edge2);
    FrontEdge3* addFrontEdge3(const Vertex3* vert0, const Vertex3* vert1);
    void addFrontEdge3(const FrontEdge3* fEdge);

    bool vertInsideFrontCheck(const Vec3& v);
    bool lineSegmentGlobalIntersectionCheck(const Vec3& v0, const Vec3& v1);
    bool segmentFrontIntersectionCheck(const Vec3& v0, const Vec3& v1);
    bool edgeGlobalIntersectionCheck(const Edge3* edge);
    bool edgeIntersectionCheck(FrontEdge3* fEdge);
    bool facetIntersectionCheck(const Vertex3* v0, const Vertex3* v1, const Vec3& v2);
    bool facetIntersectionCheck(const Vertex3* v0, const Vertex3* v1, const Vertex3* v2);
    bool facetsIntersectionCheck(FrontEdge3* fEdge);
    bool insideSimplex3Check(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3, const Vec3& vert);
    bool anyVertInsidePotentialSimp3Check(FrontEdge3* fEdge);
    bool frontSplitCheck(FrontEdge3* fEdge);
    bool parallelFacetsCheck(FrontEdge3* fEdge) const;

    static double computeMinEdgesLength(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3);
    static double computeMaxEdgesLength(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3);
    static std::pair<double, double> computeMinMaxEdgesLengths(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3);
    static std::pair<double, double> computeMinMaxEdgesSqrLengths(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3);
    static double computeSimp3SimpleQuality(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3);
    static double computeSimp3SimpleSqrQuality(const Vec3& p0, const Vec3& p1, const Vec3& p2, const Vec3& p3);

    //FrontEdge3* currentFrontEdge(double maxExCos) const;
    FrontEdge3* currentFrontEdge(double maxCompl) const;
    bool exhaustWithoutNewVertPriorityPredicate(FrontEdge3* fEdge);
    bool exhaustWithNewVertPriorityPredicate(FrontEdge3* fEdge);
    ExhaustType computeExhaustionTypeQualityPriority(
        FrontEdge3* fEdge,
        FrontFacet3** out_withNWFrontFacet = nullptr, Vec3** out_withNWNewVertPos = nullptr);

    Vec3 computeNormalInSimp3(FrontFacet3* fFacet, const Vec3& oppVertPos);
    Vec3 computeNormalInSimp3(FrontFacet3* fFacet, Edge3* oneOfRemainingEdges);

    void exhaustWithoutNewVertOppEdgeExists(FrontEdge3* fEdge, FrontEdge3* oppEdge);
    void exhaustWithoutNewVertOppEdgeDontExists(FrontEdge3* fEdge);
    void exhaustWithoutNewVert(FrontEdge3* fEdge, bool doesOppEdgeExists = true, FrontEdge3* oppEdge = nullptr);


    bool tryComputeNewVertPosType3(FrontFacet3* fFacet, Vec3& out_pos);
    bool tryComputeNewVertPosType2(
        int smallAngleIndex0, double angleCos0,
        int smallAngleIndex1, double angleCos1,
        FrontFacet3* frontFacet, Vec3& out_pos);
    bool tryComputeNewVertPosType1(
        int smallAngleIndex, double angleCos,
        FrontFacet3* fFacet, Vec3& out_pos);
    bool tryComputeNewVertPosType0(FrontFacet3* fFacet, Vec3& out_pos);
    bool tryComputeNewVertPos(FrontFacet3* fFacet, Vec3& out_pos);

    FrontFacet3* chooseFacetForExhaustionWithNewVert(FrontEdge3* fEdge);
    void         exhaustWithNewVert(FrontFacet3* fFacet, Vec3 vertPos);

    bool tryExhaustWithoutNewVert(FrontEdge3* fEdge, bool doesOppEdgeExists = true, FrontEdge3* oppEdge = nullptr);
    bool tryExhaustWithNewVert   (FrontEdge3* fEdge);

    bool globalIntersectionCheck();
    bool shellContainsVert(const Vertex3* vert);

    void processLastSimp3();
    void processAngles();

    void smoothMesh(int nIterations);
    void smoothNotFinisedMesh(int nIterations);
    void smoothFront(int nIterations);
    void smoothAroundFrontVert(Vertex3* frontVert);
};