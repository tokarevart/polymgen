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

void Edge3::Flip(
	vector<unique_ptr<Edge3>*> &edges, 
	vector<unique_ptr<Facet3>*> &facets)
{
	unique_ptr<Vertex3>* around_nodes[2];
	unique_ptr<Facet3>* around_facets[2];
	Find2FacetsAround(facets, around_facets[0], around_facets[1]);
	around_nodes[0] = (*around_facets[0])->FindVertexNotIncludedInEdge(*this);
	around_nodes[1] = (*around_facets[1])->FindVertexNotIncludedInEdge(*this);
	unique_ptr<Edge3>* new_edge = (new Edge3(**around_nodes[0], **around_nodes[1]))->GetPtrToUniquePtr();
	edges.push_back(new_edge);

	facets.push_back(
		(new Facet3(
			**(*around_facets[0])->FindEdge(**around_nodes[0], **vertexes[0]),
			**(*around_facets[1])->FindEdge(**around_nodes[1], **vertexes[0]),
			**new_edge))
		->GetPtrToUniquePtr());
	facets.push_back(
		(new Facet3(
			**(*around_facets[0])->FindEdge(**around_nodes[0], **vertexes[1]),
			**(*around_facets[1])->FindEdge(**around_nodes[1], **vertexes[1]),
			**new_edge))
		->GetPtrToUniquePtr());

	delete around_facets[0]->release();
	delete around_facets[1]->release();

	delete GetPtrToUniquePtr()->release();
}

const bool Edge3::FlipIfNeeded(
	vector<unique_ptr<Edge3>*> &edges, 
	vector<unique_ptr<Facet3>*> &facets)
{
	if (NeedToFlip(facets))
		Flip(edges, facets);

	return false;
}

void Edge3::Find2FacetsAround(
	const vector<unique_ptr<Facet3>*> &facets, 
	unique_ptr<Facet3>* &facet0, 
	unique_ptr<Facet3>* &facet1)
{
	bool not_found_yet = true;
	for (auto facet : facets)
	{
		if (*facet && 
			(*facet)->IsContaining(*this))
		{
			if (not_found_yet)
			{
				facet0 = facet;
				not_found_yet = false;
			}
			else
			{
				facet1 = facet;
				return;
			}
		}
	}
}

unique_ptr<Vertex3>* Edge3::FindFacetVertexNotBelongsToEdge(Facet3& facet, Edge3& edge)
{
	for (auto &facet_edge : facet.edges)
		if (facet_edge->get() != &edge)
			for (auto &vertex : (*facet_edge)->vertexes)
				if (!edge.IsContaining(**vertex))
					return vertex;

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

void Edge3::MakeTwoInstead(
	vector<unique_ptr<Facet3>*> &facets, 
	vector<unique_ptr<Edge3>*> &edges, 
	vector<unique_ptr<Vertex3>*> &verts)//, vector<ShellEdge3*>& shellEdges, vector<ShellFacet3*>& shellFacets)
{
	if (!*vertexes[0] || !*vertexes[1])
		throw std::exception("Something went wrong in Edge3::MakeTwoInstead");

	unique_ptr<Vertex3>* inner_vert = (new Vertex3((*vertexes[0])->GetPosition() + 0.5 * (**vertexes[1] - **vertexes[0])))->GetPtrToUniquePtr();
	verts.push_back(inner_vert);

	unique_ptr<Facet3>* old_facets[2];
	Find2FacetsAround(facets, old_facets[0], old_facets[1]);

	unique_ptr<Edge3>* edge_halfs[2];
	edge_halfs[0] = (new Edge3(**vertexes[0], **inner_vert))->GetPtrToUniquePtr();
	edge_halfs[1] = (new Edge3(**inner_vert, **vertexes[1]))->GetPtrToUniquePtr();

	unique_ptr<Vertex3>* not_belongs_to_edge_vert;
	unique_ptr<Edge3>* dividing_edge;
	for (auto &facet : old_facets)
	{
		//not_belongs_to_edge_vert = FindFacetVertexNotBelongsToEdge(**facet, *this);
		not_belongs_to_edge_vert = (*facet)->FindVertexNotIncludedInEdge(*this);
		dividing_edge = (new Edge3(**inner_vert, **not_belongs_to_edge_vert))->GetPtrToUniquePtr();
		edges.push_back(dividing_edge);

		facets.push_back(
			(new Facet3(
				**(*facet)->FindEdge(**vertexes[0], **not_belongs_to_edge_vert),
				**edge_halfs[0],
				**dividing_edge))
			->GetPtrToUniquePtr());

		facets.push_back(
			(new Facet3(
				**(*facet)->FindEdge(**vertexes[1], **not_belongs_to_edge_vert),
				**edge_halfs[1],
				**dividing_edge))
			->GetPtrToUniquePtr());

		delete facet->release();
	}

	edges.push_back(edge_halfs[0]);
	edges.push_back(edge_halfs[1]);

	delete GetPtrToUniquePtr()->release();
}

const bool Edge3::IsContaining(const Vertex3& vertex) const
{
	if (vertexes[0]->get() == &vertex ||
		vertexes[1]->get() == &vertex)
		return true;

	return false;
}

const bool Edge3::BelongsToShell()
{
	if (((*vertexes[0])->belongsToShellVertex ||
		 (*vertexes[0])->belongsToShellEdge   ||
		 (*vertexes[0])->belongsToShellFacet) &&
		((*vertexes[1])->belongsToShellVertex ||
		 (*vertexes[1])->belongsToShellEdge   ||
		 (*vertexes[1])->belongsToShellFacet))
		return true;

	return false;
}

const bool Edge3::NeedToFlip(const vector<unique_ptr<Facet3>*> &facets)
{
	unique_ptr<Vertex3>* around_nodes[2];
	unique_ptr<Facet3>* around_facets[2];
	Find2FacetsAround(
		facets, 
		around_facets[0], 
		around_facets[1]);
	around_nodes[0] = (*around_facets[0])->FindVertexNotIncludedInEdge(*this);
	around_nodes[1] = (*around_facets[1])->FindVertexNotIncludedInEdge(*this);

	double alpha = acos(Vector3::Cos(
		**vertexes[0] - **around_nodes[0], 
		**vertexes[1] - **around_nodes[0]));
	double beta = acos(Vector3::Cos(
		**vertexes[0] - **around_nodes[1], 
		**vertexes[1] - **around_nodes[1]));

	if (alpha + beta > PI)
		return true;

	return false;
}

//void Edge3::DestroyIfNoLinks()
//{
//	if (inclInFacets.empty())
//		delete _uniquePtr->release();
//}

Edge3::Edge3() : unique_ptr_helper<Edge3>(this) 
{
	vertexes[0] = nullptr;
	vertexes[1] = nullptr;
}

Edge3::Edge3(Vertex3& vertex0, Vertex3& vertex1) : unique_ptr_helper<Edge3>(this)
{
	vertexes[0] = vertex0.GetPtrToUniquePtr();
	vertexes[1] = vertex1.GetPtrToUniquePtr();
}

Edge3::~Edge3() {}

double FrontEdge3::Angle(const vector<unique_ptr<FrontFacet3>*> &frontFacets)
{
	unique_ptr<FrontFacet3>* f_facets[2];
	FindFrontFacetsAround(frontFacets, f_facets[0], f_facets[1]);

	Vector3 edge_inner_vec = (*edge->vertexes[0])->GetPosition() + (**edge->vertexes[0] - **edge->vertexes[1]);
	Vector3 facets_inner_vecs[2];
	facets_inner_vecs[0] = (*(*f_facets[0])->facet->FindVertexNotIncludedInEdge(*edge))->GetPosition() - edge_inner_vec;
	facets_inner_vecs[1] = (*(*f_facets[1])->facet->FindVertexNotIncludedInEdge(*edge))->GetPosition() - edge_inner_vec;
	Vector3 shell_inside_test = facets_inner_vecs[0] + facets_inner_vecs[1];

	int intersects_num = 0;
	for (auto &f_facet : frontFacets)
	{
		if (!*f_facet ||
			f_facet == f_facets[0] ||
			f_facet == f_facets[1])
			continue;

		if (Vector3::RayIntersectTriangle(
				edge_inner_vec,
				shell_inside_test,
				(*(*(*f_facets[0])->facet->edges[0])->vertexes[0])->GetPosition(),
				(*(*(*f_facets[0])->facet->edges[0])->vertexes[0])->GetPosition(),
				(*(*f_facets[0])->facet->FindVertexNotIncludedInEdge(**(*f_facets[0])->facet->edges[0]))->GetPosition()))
			intersects_num++;
	}

	double min_angle = Vector3::Cos(facets_inner_vecs[0], facets_inner_vecs[1]);
	
	return intersects_num % 2 == 1 ?
		min_angle :
		PI + PI - min_angle;
}

void FrontEdge3::FindFrontFacetsAround(
	const vector<unique_ptr<FrontFacet3>*> &frontFacets, 
	unique_ptr<FrontFacet3>* &facet0, 
	unique_ptr<FrontFacet3>* &facet1)
{
	bool not_found_yet = true;
	for (auto f_facet : frontFacets)
	{
		if (*f_facet &&
			(*f_facet)->facet->IsContaining(*edge))
		{
			if (not_found_yet)
			{
				facet0 = f_facet;
				not_found_yet = false;
			}
			else
			{
				facet1 = f_facet;
				return;
			}
		}
	}
}

FrontEdge3::FrontEdge3(Edge3 &edge) : edge(&edge), unique_ptr_helper<FrontEdge3>(this) {}

FrontEdge3::~FrontEdge3() {}