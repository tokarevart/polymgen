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
		//const vector<unique_ptr<Vertex3>*> &vertexes, 
		const vector<unique_ptr<Edge3>*> &edges,
		const vector<unique_ptr<Facet3>*> &facets);

	unique_ptr<FrontFacet3>* FindFrontFacet(unique_ptr<Facet3>* &facet);
	unique_ptr<FrontEdge3>* FindFrontEdge(const unique_ptr<Vertex3>* &v0, const unique_ptr<Vertex3>* &v1);
	unique_ptr<FrontEdge3>* FindFrontEdge(unique_ptr<Edge3>* &edge);

	// Adds new front facet and corresponding facet.
	void AddFrontFacet(
		unique_ptr<Edge3>* &edge0, 
		unique_ptr<Edge3>* &edge1, 
		unique_ptr<Edge3>* &edge2);

	const bool LineSegmentIntersectionCheck(const Vector3& v0, const Vector3& v1);
	const bool EdgeIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool FrontSplitCheck(const unique_ptr<FrontEdge3>* frontEdge);
	const bool FrontContainsOfOnly1Simplex3();

	void CreateSimplex3WithoutNewVertex(unique_ptr<FrontEdge3>* frontEdge);
	void CreateSimplex3WithNewVertex(unique_ptr<FrontEdge3>* frontEdge);

	void ProcessSmallAngles();
	void ProcessMediumAngles();
	void ProcessLargeAngles();
	void TriangulateVolume(const double preferredLength);

	Crystallite3();
	~Crystallite3();
};