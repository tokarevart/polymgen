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

	const bool EdgeIntersectionCheck(const Vector3& v0, const Vector3& v1);
	const bool FrontSplitCheck(const unique_ptr<FrontEdge3>* frontEdge);

	void ProcessSmallAngles();
	void ProcessMediumAngles();
	void ProcessLargeAngles();
	void TriangulateVolume(const double preferredLength);

	// Calculate based on position in space.
	//const bool Contains(const Vertex3& vertex) const;

	Crystallite3();
	~Crystallite3();
};