#include "Facet3.h"


const bool Facet3::IntersectsBy(
	const Vector3 &origin, 
	const Vector3 &dir)
{
	return Vector3::RayIntersectTriangle(
		origin, 
		dir, 
		(*(*edges[0])->vertexes[0])->GetPosition(), 
		(*(*edges[0])->vertexes[1])->GetPosition(), 
		(*FindVertexNotIncludedInEdge(**edges[0]))->GetPosition());
}

unique_ptr<Vertex3>* Facet3::FindVertexNotIncludedInEdge(const Edge3 &edge) const
{
	for (auto &facet_edge : edges)
	{
		if (&edge != facet_edge->get())
		{
			if (!edge.IsContaining(**(*facet_edge)->vertexes[0]))
				return (*facet_edge)->vertexes[0];
			else
				return (*facet_edge)->vertexes[1];
		}
	}

	return nullptr;
}

unique_ptr<Edge3>* Facet3::FindEdgeNotContainingVertex(const Vertex3 &vert) const
{
	for (auto &edge : edges)
		if (!(*edge)->IsContaining(vert))
			return edge;

	return nullptr;
}

unique_ptr<Edge3>* Facet3::FindEdge(const Vertex3 &vert0, const Vertex3 &vert1)
{
	for (auto &edge : edges)
		if (((*edge)->vertexes[0]->get() == &vert0 &&
			 (*edge)->vertexes[1]->get() == &vert1) ||
			((*edge)->vertexes[0]->get() == &vert1 &&
			 (*edge)->vertexes[1]->get() == &vert0))
			return edge;
	
	return nullptr;
}

unique_ptr<Edge3>* Facet3::MinEdge()
{
	if ((*edges[0])->SqrMagnitude() < (*edges[1])->SqrMagnitude())
	{
		if ((*edges[0])->SqrMagnitude() < (*edges[2])->SqrMagnitude())
			return edges[0];
		else
			return edges[2];
	}
	else
	{
		if ((*edges[1])->SqrMagnitude() < (*edges[2])->SqrMagnitude())
			return edges[1];
		else
			return edges[2];
	}

	return nullptr;
}

unique_ptr<Edge3>* Facet3::MaxEdge()
{
	if ((*edges[0])->SqrMagnitude() > (*edges[1])->SqrMagnitude())
	{
		if ((*edges[0])->SqrMagnitude() > (*edges[2])->SqrMagnitude())
			return edges[0];
		else
			return edges[2];
	}
	else
	{
		if ((*edges[1])->SqrMagnitude() > (*edges[2])->SqrMagnitude())
			return edges[1];
		else
			return edges[2];
	}

	return nullptr;
}

const bool Facet3::IsContaining(const Edge3 &edge) const
{
	for (auto &edge_ : edges)
		if (edge_->get() == &edge)
			return true;

	return false;
}

const bool Facet3::IsContaining(const Vertex3 &vert) const
{
	for (auto &edge : edges)
		if ((*edge)->IsContaining(vert))
			return true;

	return false;
}

Facet3::Facet3() : unique_ptr_helper<Facet3>(this) 
{
	for (auto &edge : edges)
		edge = nullptr;
}

Facet3::Facet3(Edge3 &edge0, Edge3 &edge1, Edge3 &edge2) : unique_ptr_helper<Facet3>(this)
{
	edges[0] = edge0.GetPtrToUniquePtr();
	edges[1] = edge1.GetPtrToUniquePtr();
	edges[2] = edge2.GetPtrToUniquePtr();
}

Facet3::Facet3(Vertex3 &vert0, Vertex3 &vert1, Vertex3 &vert2) : unique_ptr_helper<Facet3>(this)
{
	edges[0] = (new Edge3(vert0, vert1))->GetPtrToUniquePtr();
	edges[1] = (new Edge3(vert1, vert2))->GetPtrToUniquePtr();
	edges[2] = (new Edge3(vert2, vert0))->GetPtrToUniquePtr();
}

Facet3::~Facet3() {}

FrontFacet3::FrontFacet3(Facet3* facet) : facet(facet), unique_ptr_helper<FrontFacet3>(this) {}

FrontFacet3::~FrontFacet3() {}