#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;

class Facet3 : public unique_ptr_helper<Facet3>
{
public:
	unique_ptr<Edge3>* edges[3];

	unique_ptr<Vertex3>* FindVertexNotIncludedInEdge(const Edge3& edge) const;
	unique_ptr<Edge3>* FindEdge(const Vertex3& vertex0, const Vertex3& vertex1);
	unique_ptr<Edge3>* MinEdge();
	unique_ptr<Edge3>* MaxEdge();

	const bool IsContaining(const Edge3& edge) const;
	const bool IsContaining(const Vertex3& vertex) const;

	Facet3();
	Facet3(Edge3& edge0, Edge3& edge1, Edge3& edge2);
	Facet3(Vertex3& vertex0, Vertex3& vertex1, Vertex3& vertex2);
	~Facet3();
};

class FrontFacet3 : public unique_ptr_helper<FrontFacet3>
{
public:
	Facet3* facet;

	const Vector3 Normal(vector<unique_ptr<FrontFacet3>*> &front_facets);

	FrontFacet3(Facet3 &facet);
	~FrontFacet3();
};