#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

using std::vector;

class Crystallite3
{
public:
	vector<ShellFacet3*> shellFacets;
	vector<ShellEdge3*> shellEdges;

	//vector<unique_ptr<Vertex3>*> crysInnerVerts;
	//vector<unique_ptr<Facet3>*> innerFacets;

	vector<unique_ptr<FrontFacet3>*> frontFacets;
	vector<unique_ptr<FrontEdge3>*> frontEdges;

	const bool ShellContainsVertex(const unique_ptr<Vertex3>* &vert);
	void SetStartFront(vector<unique_ptr<Vertex3>*> &free_vertexes, vector<unique_ptr<Edge3>*> &free_edges, vector<unique_ptr<Facet3>*> &free_facets);

	// Calculate based on position in space.
	const bool Contains(const Vertex3& vertex) const;

	Crystallite3();
	~Crystallite3();
};