#include "Edge3.h"
#include <algorithm>
#include <iostream>

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

void Edge3::MakeTwoInstead(
	vector<unique_ptr<Facet3>*> &facets, 
	vector<unique_ptr<Edge3>*> &edges, 
	vector<unique_ptr<Vertex3>*> &verts)//, vector<ShellEdge3*>& shellEdges, vector<ShellFacet3*>& shellFacets)
{
	if (!*vertexes[0] || !*vertexes[1])
		throw std::logic_error("Something went wrong in Edge3::MakeTwoInstead");

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

double FrontEdge3::AngleExCos(const vector<unique_ptr<FrontFacet3>*>& frontFacets)
{
	if (!needProcessing)
		return exCos;

	unique_ptr<FrontFacet3>* f_facets[2];
	FindFrontFacetsAround(frontFacets, f_facets[0], f_facets[1]);

	Vector3 edge_inner_vec = 0.5 * ((*edge->vertexes[0])->GetPosition() + (*edge->vertexes[1])->GetPosition());
	Vector3 opp_verts_poses[2];
	opp_verts_poses[0] = (*(*f_facets[0])->facet->FindVertexNotIncludedInEdge(*edge))->GetPosition();
	opp_verts_poses[1] = (*(*f_facets[1])->facet->FindVertexNotIncludedInEdge(*edge))->GetPosition();
	Vector3 facets_inner_vecs[2];
	facets_inner_vecs[0] = opp_verts_poses[0] - Vector3::Project(opp_verts_poses[0], (*edge->vertexes[0])->GetPosition(), (*edge->vertexes[1])->GetPosition());
	facets_inner_vecs[1] = opp_verts_poses[1] - Vector3::Project(opp_verts_poses[1], (*edge->vertexes[0])->GetPosition(), (*edge->vertexes[1])->GetPosition());

	double in_the_same_plane_dot = Vector3::MixedProduct(facets_inner_vecs[0], facets_inner_vecs[1], **edge->vertexes[1] - **edge->vertexes[0]);
	if (in_the_same_plane_dot < 1e-6 && in_the_same_plane_dot > -1e-6)
		return -1.0;

	Vector3 shell_inside_test = facets_inner_vecs[0] + facets_inner_vecs[1];

	Vector3 add_for_correct_intersect = Vector3(2.1632737147, 1.488313178, -0.71123534278) * shell_inside_test.Magnitude() * 1e-3; //Vector3::CrossProduct(opp_verts_poses[0] - edge_inner_vec, opp_verts_poses[1] - edge_inner_vec).Normalize() * shell_inside_test.Magnitude() * 1e-3;
	shell_inside_test += add_for_correct_intersect;

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
			(*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition(),
			(*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition(),
			(*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition()))
			intersects_num++;
	}

	needProcessing = false;
	return exCos = intersects_num % 2 == 1 ?
		Vector3::Cos(facets_inner_vecs[0], facets_inner_vecs[1]) :
		-2.0 - Vector3::Cos(facets_inner_vecs[0], facets_inner_vecs[1]);
}

double FrontEdge3::AngleCos(bool &out_isConcave, const vector<unique_ptr<FrontFacet3>*> &frontFacets)
{
	unique_ptr<FrontFacet3>* f_facets[2];
	FindFrontFacetsAround(frontFacets, f_facets[0], f_facets[1]);

	Vector3 edge_inner_vec = 0.5 * ((*edge->vertexes[0])->GetPosition() + (*edge->vertexes[1])->GetPosition());
	Vector3 opp_verts_poses[2];
	opp_verts_poses[0] = (*(*f_facets[0])->facet->FindVertexNotIncludedInEdge(*edge))->GetPosition();
	opp_verts_poses[1] = (*(*f_facets[1])->facet->FindVertexNotIncludedInEdge(*edge))->GetPosition();
	Vector3 facets_inner_vecs[2];
	facets_inner_vecs[0] = opp_verts_poses[0] - Vector3::Project(opp_verts_poses[0], (*edge->vertexes[0])->GetPosition(), (*edge->vertexes[1])->GetPosition());
	facets_inner_vecs[1] = opp_verts_poses[1] - Vector3::Project(opp_verts_poses[1], (*edge->vertexes[0])->GetPosition(), (*edge->vertexes[1])->GetPosition());
	
	double in_the_same_plane_dot = Vector3::DotProduct(Vector3::CrossProduct(facets_inner_vecs[0], facets_inner_vecs[1]), **edge->vertexes[1] - **edge->vertexes[0]);
	if (in_the_same_plane_dot < 1e-6 && in_the_same_plane_dot > -1e-6)
		return -1.0;

	Vector3 shell_inside_test = facets_inner_vecs[0] + facets_inner_vecs[1];

	Vector3 add_for_correct_intersect = Vector3(2.1632737147, 1.488313178, -0.71123534278) * shell_inside_test.Magnitude() * 1e-3; //Vector3::CrossProduct(opp_verts_poses[0] - edge_inner_vec, opp_verts_poses[1] - edge_inner_vec).Normalize() * shell_inside_test.Magnitude() * 1e-3;
	shell_inside_test += add_for_correct_intersect;

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
				(*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition(),
				(*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition(),
				(*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition()))
			intersects_num++;
	}

	out_isConcave = intersects_num % 2 == 1;
	return Vector3::Cos(facets_inner_vecs[0], facets_inner_vecs[1]);
}

double FrontEdge3::Angle(const vector<unique_ptr<FrontFacet3>*> &frontFacets)
{
	unique_ptr<FrontFacet3>* f_facets[2];
	FindFrontFacetsAround(frontFacets, f_facets[0], f_facets[1]);

	Vector3 edge_inner_vec = 0.5 * ((*edge->vertexes[0])->GetPosition() + (*edge->vertexes[1])->GetPosition());
	Vector3 opp_verts_poses[2];
	opp_verts_poses[0] = (*(*f_facets[0])->facet->FindVertexNotIncludedInEdge(*edge))->GetPosition();
	opp_verts_poses[1] = (*(*f_facets[1])->facet->FindVertexNotIncludedInEdge(*edge))->GetPosition();
	Vector3 facets_inner_vecs[2];
	facets_inner_vecs[0] = opp_verts_poses[0] - Vector3::Project(opp_verts_poses[0], (*edge->vertexes[0])->GetPosition(), (*edge->vertexes[1])->GetPosition());
	facets_inner_vecs[1] = opp_verts_poses[1] - Vector3::Project(opp_verts_poses[1], (*edge->vertexes[0])->GetPosition(), (*edge->vertexes[1])->GetPosition());
	double in_the_same_plane_dot = Vector3::DotProduct(Vector3::CrossProduct(facets_inner_vecs[0], facets_inner_vecs[1]), **edge->vertexes[1] - **edge->vertexes[0]);
	if (in_the_same_plane_dot < 1e-6 && in_the_same_plane_dot > -1e-6)
		return PI;
	Vector3 shell_inside_test = facets_inner_vecs[0] + facets_inner_vecs[1];

	Vector3 add_for_correct_intersect = Vector3(2.1632737147, 1.488313178, -0.71123534278) * shell_inside_test.Magnitude() * 1e-3; //Vector3::CrossProduct(opp_verts_poses[0] - edge_inner_vec, opp_verts_poses[1] - edge_inner_vec).Normalize() * shell_inside_test.Magnitude() * 1e-3;
	shell_inside_test += add_for_correct_intersect;

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
				(*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition(),
				(*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition(),
				(*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition()))
			intersects_num++;
	}

	double min_angle = acos(Vector3::Cos(facets_inner_vecs[0], facets_inner_vecs[1]));
	/*if (intersects_num % 2 == 0)
		std::cout << "wrong intersect";*/
	return intersects_num % 2 == 1 ?
		min_angle :
		PI + PI - min_angle;
}

void FrontEdge3::FindFrontFacetsAround(
	const vector<unique_ptr<FrontFacet3>*> &frontFacets, 
	unique_ptr<FrontFacet3>* &out_frontFacet0,
	unique_ptr<FrontFacet3>* &out_frontFacet1)
{
	bool not_found_yet = true;
	for (auto f_facet : frontFacets)
	{
		if (*f_facet &&
			(*f_facet)->facet->IsContaining(*edge))
		{
			if (not_found_yet)
			{
				out_frontFacet0 = f_facet;
				not_found_yet = false;
			}
			else
			{
				out_frontFacet1 = f_facet;
				return;
			}
		}
	}

	throw std::logic_error("Something went wrong in FrontEdge3::FindFrontFacetsAround");
}

void FrontEdge3::FindOppositeVertexes(
	const vector<unique_ptr<FrontFacet3>*> &frontFacets, 
	const vector<unique_ptr<FrontEdge3>*> &frontEdges, 
	unique_ptr<Vertex3>* &out_vert0,
	unique_ptr<Vertex3>* &out_vert1)
{
	unique_ptr<FrontFacet3>* around_f_facets[2];
	FindFrontFacetsAround(
		frontFacets, 
		around_f_facets[0], 
		around_f_facets[1]);

	out_vert0 = (*around_f_facets[0])->facet->FindVertexNotIncludedInEdge(*edge);
	out_vert1 = (*around_f_facets[1])->facet->FindVertexNotIncludedInEdge(*edge);
}

unique_ptr<FrontEdge3>* FrontEdge3::FindOppositeFrontEdge(
	const vector<unique_ptr<FrontFacet3>*> &frontFacets, 
	const vector<unique_ptr<FrontEdge3>*> &frontEdges)
{
	unique_ptr<Vertex3>* opp_verts[2];
	FindOppositeVertexes(
		frontFacets,
		frontEdges,
		opp_verts[0],
		opp_verts[1]);

	for (auto f_edge : frontEdges)
	{
		if (*f_edge &&
			(((*f_edge)->edge->vertexes[0] == opp_verts[0] &&
			  (*f_edge)->edge->vertexes[1] == opp_verts[1]) ||
			 ((*f_edge)->edge->vertexes[1] == opp_verts[0] &&
			  (*f_edge)->edge->vertexes[0] == opp_verts[1])))
			return f_edge;
	}

	return nullptr;
}

FrontEdge3::FrontEdge3(Edge3* edge) : edge(edge), unique_ptr_helper<FrontEdge3>(this) {}

FrontEdge3::~FrontEdge3() {}