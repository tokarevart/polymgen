#include "Facet3.h"


unique_ptr<Vertex3>* Facet3::FindVertexNotIncludedInEdge(const Edge3& edge) const
{
	for (auto &facet_edge : edges)
	{
		if (&edge != facet_edge->get())
		{
			if (!edge.IsContaining(**(*facet_edge)->vertexes[0]))
			{
				return (*facet_edge)->vertexes[0];
			}
			else
			{
				return (*facet_edge)->vertexes[1];
			}
		}
	}

	return nullptr;
}

unique_ptr<Edge3>* Facet3::FindEdge(const Vertex3& vertex0, const Vertex3& vertex1)
{
	for (auto &edge : edges)
	{
		if (((*edge)->vertexes[0]->get() == &vertex0 &&
			 (*edge)->vertexes[1]->get() == &vertex1) ||
			((*edge)->vertexes[0]->get() == &vertex1 &&
			 (*edge)->vertexes[1]->get() == &vertex0))
		{
			return edge;
		}
	}
	
	return nullptr;
}

unique_ptr<Edge3>* Facet3::MinEdge()
{
	if ((*edges[0])->SqrMagnitude() < (*edges[1])->SqrMagnitude())
	{
		if ((*edges[0])->SqrMagnitude() < (*edges[2])->SqrMagnitude())
		{
			return edges[0];
		}
		else
		{
			return edges[2];
		}
	}
	else
	{
		if ((*edges[1])->SqrMagnitude() < (*edges[2])->SqrMagnitude())
		{
			return edges[1];
		}
		else
		{
			return edges[2];
		}
	}

	return nullptr;
}

unique_ptr<Edge3>* Facet3::MaxEdge()
{
	if ((*edges[0])->SqrMagnitude() > (*edges[1])->SqrMagnitude())
	{
		if ((*edges[0])->SqrMagnitude() > (*edges[2])->SqrMagnitude())
		{
			return edges[0];
		}
		else
		{
			return edges[2];
		}
	}
	else
	{
		if ((*edges[1])->SqrMagnitude() > (*edges[2])->SqrMagnitude())
		{
			return edges[1];
		}
		else
		{
			return edges[2];
		}
	}

	return nullptr;
}

const bool Facet3::IsContaining(const Edge3& edge) const
{
	for (int i = 0; i < 3; i++)
	{
		if (edges[i]->get() == &edge)
		{
			return true;
		}
	}

	return false;
}

const bool Facet3::IsContaining(const Vertex3& vertex) const
{
	for (auto &edge : edges)
	{
		if ((*edge)->IsContaining(vertex))
		{
			return true;
		}
	}

	return false;
}

Facet3::Facet3() : unique_ptr_helper<Facet3>(this) 
{
	for (auto &edge : edges)
	{
		edge = nullptr;
	}
}

Facet3::Facet3(Edge3& edge0, Edge3& edge1, Edge3& edge2) : unique_ptr_helper<Facet3>(this)
{
	edges[0] = edge0.GetPtrToUniquePtr();
	edges[1] = edge1.GetPtrToUniquePtr();
	edges[2] = edge2.GetPtrToUniquePtr();

	(*edges[0])->inclInFacets.push_back(GetPtrToUniquePtr());
	(*edges[1])->inclInFacets.push_back(GetPtrToUniquePtr());
	(*edges[2])->inclInFacets.push_back(GetPtrToUniquePtr());

	for (auto &edge : edges)
	{
		for (auto &vertex : (*edge)->vertexes)
		{
			if (std::find(
					(*vertex)->inclInFacets.begin(),
					(*vertex)->inclInFacets.end(), 
					GetPtrToUniquePtr())
				== (*vertex)->inclInFacets.end())
			{
				(*vertex)->inclInFacets.push_back(GetPtrToUniquePtr());
			}
		}
	}
}

Facet3::Facet3(Vertex3& vertex0, Vertex3& vertex1, Vertex3& vertex2) : unique_ptr_helper<Facet3>(this)
{
	edges[0] = (new Edge3(vertex0, vertex1))->GetPtrToUniquePtr();
	edges[1] = (new Edge3(vertex1, vertex2))->GetPtrToUniquePtr();
	edges[2] = (new Edge3(vertex2, vertex0))->GetPtrToUniquePtr();

	(*edges[0])->inclInFacets.push_back(GetPtrToUniquePtr());
	(*edges[1])->inclInFacets.push_back(GetPtrToUniquePtr());
	(*edges[2])->inclInFacets.push_back(GetPtrToUniquePtr());

	for (auto &edge : edges)
	{
		for (auto &vertex : (*edge)->vertexes)
		{
			if (std::find(
					(*vertex)->inclInFacets.begin(),
					(*vertex)->inclInFacets.end(),
					GetPtrToUniquePtr())
				== (*vertex)->inclInFacets.end())
			{
				(*vertex)->inclInFacets.push_back(GetPtrToUniquePtr());
			}
		}
	}
}

Facet3::~Facet3()
{
	for (auto &edge : edges)
	{
		if (*edge)
		{
			for (auto &vertex : (*edge)->vertexes)
			{
				if (*vertex)
				{
					(*vertex)->inclInFacets.remove(GetPtrToUniquePtr());
				}
			}

			(*edge)->inclInFacets.remove(GetPtrToUniquePtr());
			(*edge)->DestroyIfNoLinks();
		}
	}
}