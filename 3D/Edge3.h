#pragma once
#include <list>
#include <vector>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;
using std::list;
using std::vector;

class Edge3 : public unique_ptr_helper<Edge3>
{
public:
	unique_ptr<Vertex3>* vertexes[2];

	list<unique_ptr<Facet3>*> inclInFacets;

	const double Magnitude() const;
	const double SqrMagnitude() const;

	void Flip(vector<unique_ptr<Edge3>*>& free_edges, vector<unique_ptr<Facet3>*>& free_facets);
	const bool FlipIfNeeded(vector<unique_ptr<Edge3>*>& free_edges, vector<unique_ptr<Facet3>*>& free_facets);
	void MakeTwoInstead(vector<unique_ptr<Facet3>*>& freeFacets, vector<unique_ptr<Edge3>*>& freeEdges, vector<unique_ptr<Vertex3>*>& freeVertexes);//, vector<ShellEdge3*>& shellEdges, vector<ShellFacet3*>& shellFacets);

	const bool IsContaining(const Vertex3& vertex) const;

	const bool BelongsToShell();

	const bool NeedToFlip();


	void DestroyIfNoLinks();

	Edge3();
	Edge3(Vertex3& vertex0, Vertex3& vertex1);
	~Edge3();
};

class FrontEdge3 : public unique_ptr_helper<FrontEdge3>
{
public:
	Edge3* edge;

	void FindFrontFacetsAround(vector<unique_ptr<FrontFacet3>*> &front_facets, unique_ptr<FrontFacet3>* &facet0, unique_ptr<FrontFacet3>* &facet1);

	FrontEdge3(Edge3 &edge);
	~FrontEdge3();
};