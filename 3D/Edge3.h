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

	//list<unique_ptr<Facet3>*> inclInFacets;

	const double Magnitude() const;
	const double SqrMagnitude() const;

	void Flip(
		vector<unique_ptr<Edge3>*> &edges, 
		vector<unique_ptr<Facet3>*> &facets);
	const bool FlipIfNeeded(
		vector<unique_ptr<Edge3>*> &edges, 
		vector<unique_ptr<Facet3>*> &facets);
	void Find2FacetsAround(
		const vector<unique_ptr<Facet3>*> &facets, 
		unique_ptr<Facet3>* &facet0, 
		unique_ptr<Facet3>* &facet1);
	unique_ptr<Vertex3>* FindFacetVertexNotBelongsToEdge(
		Facet3& facet, 
		Edge3& edge);
	// Doesn't change vertexes relations.
	void MakeTwoInstead(
		vector<unique_ptr<Facet3>*> &facets, 
		vector<unique_ptr<Edge3>*> &edges, 
		vector<unique_ptr<Vertex3>*> &vertexes);//, vector<ShellEdge3*>& shellEdges, vector<ShellFacet3*>& shellFacets);

	const bool IsContaining(const Vertex3& vertex) const;

	const bool BelongsToShell();

	const bool NeedToFlip(const vector<unique_ptr<Facet3>*> &facets);

	//void DestroyIfNoLinks();

	Edge3();
	Edge3(Vertex3& vertex0, Vertex3& vertex1);
	~Edge3();
};

class FrontEdge3 : public unique_ptr_helper<FrontEdge3>
{
public:
	Edge3* edge;

	double Angle(
		const vector<unique_ptr<FrontFacet3>*> &frontFacets);
	void FindFrontFacetsAround(
		const vector<unique_ptr<FrontFacet3>*> &frontFacets, 
		unique_ptr<FrontFacet3>* &facet0, 
		unique_ptr<FrontFacet3>* &facet1);

	FrontEdge3(Edge3 &edge);
	~FrontEdge3();
};