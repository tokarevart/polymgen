#include "Edge3.h"
#include <algorithm>

#define PI 3.141592653589793

#define EPS 1e-10
#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))


const double Edge3::Magnitude() const
{
	return (**vertexes[1] - **vertexes[0]).Magnitude();
}

const double Edge3::SqrMagnitude() const
{
	return (**vertexes[1] - **vertexes[0]).SqrMagnitude();
}

void Edge3::Flip(vector<unique_ptr<Edge3>*>& free_edges, vector<unique_ptr<Facet3>*>& free_facets)
{
	unique_ptr<Vertex3>* around_nodes[2];
	around_nodes[0] = (*inclInFacets.front())->FindVertexNotIncludedInEdge(*this);
	around_nodes[1] = (*inclInFacets.back())->FindVertexNotIncludedInEdge(*this);
	unique_ptr<Edge3>* new_edge = (new Edge3(**around_nodes[0], **around_nodes[1]))->GetPtrToUniquePtr();
	free_edges.push_back(new_edge);

	free_facets.push_back(
		(new Facet3(
			**(*inclInFacets.front())->FindEdge(**around_nodes[0], **vertexes[0]),
			**(*inclInFacets.back())->FindEdge(**around_nodes[1], **vertexes[0]),
			**new_edge))
		->GetPtrToUniquePtr());
	free_facets.push_back(
		(new Facet3(
			**(*inclInFacets.front())->FindEdge(**around_nodes[0], **vertexes[1]),
			**(*inclInFacets.back())->FindEdge(**around_nodes[1], **vertexes[1]),
			**new_edge))
		->GetPtrToUniquePtr());

	delete GetPtrToUniquePtr()->release();
}

const bool Edge3::FlipIfNeeded(vector<unique_ptr<Edge3>*>& free_edges, vector<unique_ptr<Facet3>*>& free_facets)
{
	if (NeedToFlip())
	{
		Flip(free_edges, free_facets);
	}

	return false;
}

unique_ptr<Vertex3>* FindFacetVertexNotBelongsToEdge(Facet3& facet, Edge3& edge)
{
	for (auto &facet_edge : facet.edges)
	{
		if (facet_edge->get() != &edge)
		{
			for (auto &vertex : (*facet_edge)->vertexes)
			{
				if (!edge.IsContaining(**vertex))
				{
					return vertex;
				}
			}
		}
	}

	return nullptr;
}

//ShellFacet3* FindShellFacet(ShellVertex3& s_vert, ShellEdge3& s_edge, vector<ShellFacet3*>& shellFacets)
//{
//	for (auto &facet : shellFacets)
//	{
//		if (facet->IsContaining(s_vert))
//		{
//			return facet;
//		}
//	}
//
//	return nullptr;
//}
//
//ShellEdge3* FindShellEdge(ShellVertex3& s_vert0, ShellVertex3& s_vert1, vector<ShellEdge3*>& shellEdges)
//{
//	for (auto &edge : shellEdges)
//	{
//		if (edge->vertexes[0] == &s_vert0 &&
//			edge->vertexes[1] == &s_vert1 ||
//			edge->vertexes[0] == &s_vert1 &&
//			edge->vertexes[1] == &s_vert0)
//		{
//			return edge;
//		}
//	}
//
//	return nullptr;
//}

void Edge3::MakeTwoInstead(vector<unique_ptr<Facet3>*>& freeFacets, vector<unique_ptr<Edge3>*>& freeEdges, vector<unique_ptr<Vertex3>*>& freeVertexes)//, vector<ShellEdge3*>& shellEdges, vector<ShellFacet3*>& shellFacets)
{
	if (!*vertexes[0] || !*vertexes[1])
	{
		throw std::exception("Something went wrong in Edge3::MakeTwoInstead");
	}

	unique_ptr<Vertex3>* inner_vert = (new Vertex3((*vertexes[0])->GetPosition() + 0.5 * (**vertexes[1] - **vertexes[0])))->GetPtrToUniquePtr();
	//if (BelongsToShell())
	//{
	//	ShellFacet3* s_facet_buf;
	//	ShellEdge3* s_edge_buf = FindShellEdge(*(*vertexes[0])->belongsToShellVertex, *(*vertexes[1])->belongsToShellVertex, shellEdges);
	//	if ((*vertexes[0])->belongsToShellVertex &&
	//		(*vertexes[1])->belongsToShellVertex &&
	//		s_edge_buf)
	//	{
	//		(*inner_vert)->belongsToShellEdge = s_edge_buf;
	//	}
	//	else if ((*vertexes[0])->belongsToShellEdge &&
	//		(*vertexes[1])->belongsToShellVertex &&
	//		(*vertexes[0])->belongsToShellEdge->IsContaining(*(*vertexes[1])->belongsToShellVertex))
	//	{
	//		(*inner_vert)->belongsToShellEdge = (*vertexes[0])->belongsToShellEdge;
	//	}
	//	else if ((*vertexes[1])->belongsToShellEdge &&
	//		(*vertexes[0])->belongsToShellVertex &&
	//		(*vertexes[1])->belongsToShellEdge->IsContaining(*(*vertexes[0])->belongsToShellVertex))
	//	{
	//		(*inner_vert)->belongsToShellEdge = (*vertexes[1])->belongsToShellEdge;
	//	}
	//	else if ((*vertexes[0])->belongsToShellEdge == (*vertexes[1])->belongsToShellEdge &&
	//		(*vertexes[0])->belongsToShellEdge)
	//	{
	//		(*inner_vert)->belongsToShellEdge = (*vertexes[0])->belongsToShellEdge;
	//	}
	//	else if ((*vertexes[0])->belongsToShellEdge != (*vertexes[1])->belongsToShellEdge &&
	//		(*vertexes[0])->belongsToShellEdge && (*vertexes[1])->belongsToShellEdge)
	//	{
	//		if ((*vertexes[0])->belongsToShellEdge->vertexes[0] == (*vertexes[1])->belongsToShellEdge->vertexes[0] ||
	//			(*vertexes[0])->belongsToShellEdge->vertexes[0] == (*vertexes[1])->belongsToShellEdge->vertexes[1])
	//		{
	//			s_facet_buf = FindShellFacet(*(*vertexes[0])->belongsToShellEdge->vertexes[1], *(*vertexes[1])->belongsToShellEdge, shellFacets);
	//		}
	//		else
	//		{
	//			s_facet_buf = FindShellFacet(*(*vertexes[0])->belongsToShellEdge->vertexes[0], *(*vertexes[1])->belongsToShellEdge, shellFacets);
	//		}
	//		(*inner_vert)->belongsToShellFacet = s_facet_buf;
	//	}
	//	else if ((*vertexes[0])->belongsToShellEdge &&
	//		(*vertexes[1])->belongsToShellVertex &&
	//		(s_facet_buf = FindShellFacet(*(*vertexes[1])->belongsToShellVertex, *(*vertexes[0])->belongsToShellEdge, shellFacets)))
	//	{
	//		(*inner_vert)->belongsToShellFacet = s_facet_buf;
	//	}
	//	else if ((*vertexes[1])->belongsToShellEdge &&
	//		(*vertexes[0])->belongsToShellVertex &&
	//		(s_facet_buf = FindShellFacet(*(*vertexes[0])->belongsToShellVertex, *(*vertexes[1])->belongsToShellEdge, shellFacets)))
	//	{
	//		(*inner_vert)->belongsToShellFacet = s_facet_buf;
	//	}
	//	else if ((*vertexes[0])->belongsToShellVertex &&
	//		(*vertexes[1])->belongsToShellFacet)
	//	{
	//		(*inner_vert)->belongsToShellFacet = (*vertexes[1])->belongsToShellFacet;
	//	}
	//	else if ((*vertexes[1])->belongsToShellVertex &&
	//		(*vertexes[0])->belongsToShellFacet)
	//	{
	//		(*inner_vert)->belongsToShellFacet = (*vertexes[0])->belongsToShellFacet;
	//	}
	//	else if ((*vertexes[0])->belongsToShellEdge &&
	//		(*vertexes[1])->belongsToShellFacet)
	//	{
	//		(*inner_vert)->belongsToShellFacet = (*vertexes[1])->belongsToShellFacet;
	//	}
	//	else if ((*vertexes[1])->belongsToShellEdge &&
	//		(*vertexes[0])->belongsToShellFacet)
	//	{
	//		(*inner_vert)->belongsToShellFacet = (*vertexes[0])->belongsToShellFacet;
	//	}
	//	else if ((*vertexes[0])->belongsToShellFacet == (*vertexes[1])->belongsToShellFacet &&
	//		(*vertexes[0])->belongsToShellFacet)
	//	{
	//		(*inner_vert)->belongsToShellFacet = (*vertexes[0])->belongsToShellFacet;
	//	}
	//	else
	//	{
	//		throw std::exception("Something went wrong in Edge3::MakeTwoInstead");
	//	}
	//}
	freeVertexes.push_back(inner_vert);

	vector<unique_ptr<Facet3>*> oldFacets;
	for (auto &facet : inclInFacets)
	{
		if (*facet)
		{
			oldFacets.push_back(facet);
		}
	}

	unique_ptr<Edge3>* edge_halfs[2];
	edge_halfs[0] = (new Edge3(**vertexes[0], **inner_vert))->GetPtrToUniquePtr();
	edge_halfs[1] = (new Edge3(**inner_vert, **vertexes[1]))->GetPtrToUniquePtr();

	unique_ptr<Vertex3>* not_belongs_to_edge_vert;
	unique_ptr<Edge3>* dividing_edge;
	for (auto &facet : oldFacets)
	{
		not_belongs_to_edge_vert = FindFacetVertexNotBelongsToEdge(**facet, *this);
		dividing_edge = (new Edge3(**inner_vert, **not_belongs_to_edge_vert))->GetPtrToUniquePtr();
		freeEdges.push_back(dividing_edge);

		freeFacets.push_back(
			(new Facet3(
				**(*facet)->FindEdge(**vertexes[0], **not_belongs_to_edge_vert),
				**edge_halfs[0],
				**dividing_edge))
			->GetPtrToUniquePtr());

		freeFacets.push_back(
			(new Facet3(
				**(*facet)->FindEdge(**vertexes[1], **not_belongs_to_edge_vert),
				**edge_halfs[1],
				**dividing_edge))
			->GetPtrToUniquePtr());

		delete facet->release();
	}

	freeEdges.push_back(edge_halfs[0]);
	freeEdges.push_back(edge_halfs[1]);
}

const bool Edge3::IsContaining(const Vertex3& vertex) const
{
	if (vertexes[0]->get() == &vertex ||
		vertexes[1]->get() == &vertex)
	{
		return true;
	}

	return false;
}

const bool Edge3::BelongsToShell()
{
	//if (((*vertexes[0])->belongsToShellVertex && (*vertexes[1])->belongsToShellVertex) ||
	//	((*vertexes[0])->belongsToShellEdge && (*vertexes[1])->belongsToShellVertex) ||
	//	((*vertexes[1])->belongsToShellEdge && (*vertexes[0])->belongsToShellVertex) ||
	//	((*vertexes[0])->belongsToShellEdge && (*vertexes[1])->belongsToShellEdge) ||
	//	((*vertexes[0])->belongsToShellVertex && (*vertexes[1])->belongsToShellFacet) ||
	//	((*vertexes[1])->belongsToShellVertex && (*vertexes[0])->belongsToShellFacet) ||
	//	((*vertexes[0])->belongsToShellEdge && (*vertexes[1])->belongsToShellFacet) ||
	//	((*vertexes[1])->belongsToShellEdge && (*vertexes[0])->belongsToShellFacet) ||
	//	((*vertexes[0])->belongsToShellFacet  && (*vertexes[1])->belongsToShellFacet))
	if (
		((*vertexes[0])->belongsToShellVertex ||
		 (*vertexes[0])->belongsToShellEdge ||
		 (*vertexes[0])->belongsToShellFacet) &&
		((*vertexes[1])->belongsToShellVertex ||
		 (*vertexes[1])->belongsToShellEdge ||
		 (*vertexes[1])->belongsToShellFacet)
		)
	{
		return true;
	}

	return false;
}

const bool Edge3::NeedToFlip()
{
	if (inclInFacets.size() < 2)
	{
		throw std::exception("Edge included only in 1 simplex");
	}

	unique_ptr<Vertex3>* around_nodes[2];
	around_nodes[0] = (*inclInFacets.front())->FindVertexNotIncludedInEdge(*this);
	around_nodes[1] = (*inclInFacets.back())->FindVertexNotIncludedInEdge(*this);

	double alpha = acos(Vector3::Cos(**vertexes[0] - **around_nodes[0], **vertexes[1] - **around_nodes[0]));
	double beta = acos(Vector3::Cos(**vertexes[0] - **around_nodes[1], **vertexes[1] - **around_nodes[1]));

	if (alpha + beta > PI)
	{
		return true;
	}

	return false;
}

void Edge3::DestroyIfNoLinks()
{
	if (inclInFacets.empty())
	{
		delete _uniquePtr->release();
	}
}

Edge3::Edge3() : unique_ptr_helper<Edge3>(this) 
{
	for (auto &vertex : vertexes)
	{
		vertex = nullptr;
	}
}

Edge3::Edge3(Vertex3& vertex0, Vertex3& vertex1) : unique_ptr_helper<Edge3>(this)
{
	vertexes[0] = vertex0.GetPtrToUniquePtr();
	vertexes[1] = vertex1.GetPtrToUniquePtr();

	(*vertexes[0])->inclInEdges.push_back(GetPtrToUniquePtr());
	(*vertexes[0])->neighbors.push_back(vertexes[1]);

	(*vertexes[1])->inclInEdges.push_back(GetPtrToUniquePtr());
	(*vertexes[1])->neighbors.push_back(vertexes[0]);
}

Edge3::~Edge3()
{
	for (auto &nodei : vertexes)
	{
		if (*nodei)
		{
			(*nodei)->inclInEdges.remove(GetPtrToUniquePtr());
			for (auto &nodej : vertexes)
			{
				if (*nodej && (nodej != nodei))
				{
					(*nodei)->neighbors.remove(nodej);
				}
			}
			(*nodei)->DestroyIfNoLinks();
		}
	}

	for (auto &facet : inclInFacets)
	{
		if (*facet)
		{
			delete facet->release();
		}
	}
}

FrontEdge3::FrontEdge3(Edge3 &edge) : edge(&edge), unique_ptr_helper<FrontEdge3>(this) {}

FrontEdge3::~FrontEdge3() {}