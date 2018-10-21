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
	void ErasePtrsToNullptr(vector<unique_ptr<T>*> &vec);

public:
	vector<ShellFacet3*> shellFacets;
	vector<ShellEdge3*> shellEdges;

	//vector<unique_ptr<Simplex3>*> innerSimps;
	vector<unique_ptr<Facet3>*> innerFacets;
	vector<unique_ptr<Edge3>*> innerEdges;
	vector<unique_ptr<Vertex3>*> innerVerts;

	vector<unique_ptr<FrontFacet3>*> frontFacets;
	vector<unique_ptr<FrontEdge3>*> frontEdges;

	const bool ShellContainsVertex(
		const Vertex3 &vert);
	void SetStartFront(
		const vector<unique_ptr<Edge3>*> &edges,
		const vector<unique_ptr<Facet3>*> &facets);

	unique_ptr<FrontFacet3>* FindFrontFacet(unique_ptr<Facet3>* &facet);
	unique_ptr<FrontEdge3>* FindFrontEdge(const unique_ptr<Vertex3>* &v0, const unique_ptr<Vertex3>* &v1);
	unique_ptr<FrontEdge3>* FindFrontEdge(unique_ptr<Edge3>* &edge);

	// Adds new front facet and corresponding facet.
	unique_ptr<FrontFacet3>* AddFrontFacet3(
		unique_ptr<Edge3>* &edge0, 
		unique_ptr<Edge3>* &edge1, 
		unique_ptr<Edge3>* &edge2);
	unique_ptr<FrontEdge3>* AddFrontEdge3(
		unique_ptr<Vertex3>* &vert0,
		unique_ptr<Vertex3>* &vert1);
	void AddFrontEdge3(
		unique_ptr<FrontEdge3>* &frontEdge);

	const bool VertexInsideFrontCheck(const Vector3& v);
	const bool LineSegmentGlobalIntersectionCheck(const Vector3& v0, const Vector3& v1);
	const bool LineSegmentFrontIntersectionCheck(const Vector3& v0, const Vector3& v1);
	const bool EdgeGlobalIntersectionCheck(const unique_ptr<Edge3>* edge);
	const bool EdgeIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool FrontSplitCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool ParallelFacetsCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool FrontContainsOfOnly1Simplex3OrEmpty();

	void CreateSimplex3AroundEdge(unique_ptr<FrontEdge3>* frontEdge);
	Vector3 NewVertexPositionForSimplexesCreation(unique_ptr<FrontEdge3>* frontEdge, const bool calcVertDirByAngle = false);
	void Create2Simplexes3AroundEdge(unique_ptr<FrontEdge3>* frontEdge, Vector3 vertPos);

	const bool ProcessSmallAngles();
	//void ProcessSmallAngles(Polycrystalline3* polycr);
	void ProcessMediumAngles(Polycrystalline3* polycr);
	void ProcessLargeAngles();
	const bool GlobalIntersectionCheck();
	void TriangulateVolume(const double preferredLength, Polycrystalline3* polycr);

	Crystallite3();
	~Crystallite3();
};