#pragma once
#include <list>
#include <vector>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

using std::vector;
using std::unique_ptr;

class Crystallite3
{
	double _preferredLength;

	template <class T>
	void ErasePtrsToNullptr(vector<unique_ptr<T>*>& vec);

public:
	vector<ShellFacet3*> shellFacets;
	vector<ShellEdge3*> shellEdges;

	vector<unique_ptr<Simplex3>*> innerSimps;
	vector<unique_ptr<Facet3>*> innerFacets;
	vector<unique_ptr<Edge3>*> innerEdges;
	vector<unique_ptr<Vertex3>*> innerVerts;

	vector<unique_ptr<FrontFacet3>*> frontFacets;
	vector<unique_ptr<FrontEdge3>*> frontEdges;

	const bool ShellContainsVertex(
		const Vertex3& vert);
	void SetStartFront(
		const vector<unique_ptr<Edge3>*>& edges,
		const vector<unique_ptr<Facet3>*>& facets);

	ShellEdge3* FindShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const;

	unique_ptr<FrontFacet3>* FindFrontFacet(unique_ptr<Facet3>* facet);
	unique_ptr<FrontEdge3>* FindFrontEdge(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1);
	unique_ptr<FrontEdge3>* FindFrontEdge(unique_ptr<Edge3>* edge);

	// Adds new front facet and corresponding facet.
	unique_ptr<FrontFacet3>* AddFrontFacet3(
		unique_ptr<Edge3>*& edge0, 
		unique_ptr<Edge3>*& edge1, 
		unique_ptr<Edge3>*& edge2);
	unique_ptr<FrontEdge3>* AddFrontEdge3(
		unique_ptr<Vertex3>*& vert0,
		unique_ptr<Vertex3>*& vert1);
	void AddFrontEdge3(
		unique_ptr<FrontEdge3>*& frontEdge);

	const bool VertexInsideFrontCheck(const Vector3& v);
	const bool LineSegmentGlobalIntersectionCheck(const Vector3& v0, const Vector3& v1);
	const bool LineSegmentFrontIntersectionCheck(const Vector3& v0, const Vector3& v1);
	const bool EdgeGlobalIntersectionCheck(const unique_ptr<Edge3>* edge);
	const bool EdgeIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool FacetIntersectionCheck(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1, const Vector3& v2);
	const bool FacetIntersectionCheck(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1, const unique_ptr<Vertex3>* v2);
	const bool FacetsIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool InsideSimplex3Check(const Vector3& p0, const Vector3& p1, const Vector3& p2, const Vector3& p3, const Vector3& vert);
	const bool SomeVertexInsidePotentialSimplex3Check(const unique_ptr<FrontEdge3>* frontEdge);
	const bool FrontSplitCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool ParallelFacetsCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool FrontContainsOfOnly1FacetOrEmpty();

	unique_ptr<FrontEdge3>* CurrentFrontEdge(double maxExCos);
	const bool ExhaustWithoutNewVertexPriorityPredicate(unique_ptr<FrontEdge3>* currentFrontEdge);
	const bool ExhaustWithNewVertexPriorityPredicate(unique_ptr<FrontEdge3>* currentFrontEdge);
	const int  ExhaustTypeQualityPriorityCalculation(unique_ptr<FrontEdge3>* currentFrontEdge);

	void ExhaustWithoutNewVertexOppositeEdgeExists(unique_ptr<FrontEdge3>* frontEdge, unique_ptr<FrontEdge3>* oppositeEdge);
	void ExhaustWithoutNewVertexOppositeEdgeDontExists(unique_ptr<FrontEdge3>* frontEdge);
	void ExhaustWithoutNewVertex(unique_ptr<FrontEdge3>* frontEdge, const bool oppositeEdgeExistence = true, unique_ptr<FrontEdge3>* oppositeEdge = nullptr);
	const bool NewVertexPosition(unique_ptr<FrontFacet3>* frontFacet, Vector3& out_pos);
	unique_ptr<FrontFacet3>* ChooseFrontFacetForExhaustionWithNewVertex(unique_ptr<FrontEdge3>* frontEdge);
	void ExhaustWithNewVertex(unique_ptr<FrontFacet3>* frontFacet, Vector3 vertPos);
	const bool NewVertexPosition_OLD(unique_ptr<FrontEdge3>* frontEdge, Vector3& out_pos);
	void ExhaustWithNewVertex_OLD(unique_ptr<FrontEdge3>* frontEdge, Vector3 vertPos);

	bool TryExhaustWithoutNewVertex(unique_ptr<FrontEdge3>* frontEdge, const bool oppositeEdgeExistence = true, unique_ptr<FrontEdge3>* oppositeEdge = nullptr);
	bool TryExhaustWithNewVertex(unique_ptr<FrontEdge3>* frontEdge);

	const bool ProcessVerySmallAngles(Polycrystal3 *polycr);
	const bool ProcessSmallAngles(Polycrystal3 *polycr);
	const bool ProcessMediumAngles(Polycrystal3* polycr);
	void ProcessLargeAngles(Polycrystal3* polycr);
	void ProcessAngles(Polycrystal3* polycr);
	const bool GlobalIntersectionCheck();
	void SmoothMesh(int iterationsNum);
	void SmoothNotFinisedMesh(int iterationsNum);
	void SmoothFront(int iterationsNum);
	void SmoothAroundFrontVertex(unique_ptr<Vertex3>* frontVert);
	void AnalyzeMeshQuality(double& out_minQuality, double& out_avQuality);
	void TriangulateVolume(const double preferredLength, Polycrystal3* polycr);

	Crystallite3();
	~Crystallite3();
};