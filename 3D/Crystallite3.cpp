#include "Crystallite3.h"
#include <algorithm>

#define DET(a, b, c, d) \
		(a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))

#define DEG_80_IN_RADIANS 1.3962634015954636
#define DEG_160_IN_RADIANS 2.7925268031909273

const bool Crystallite3::ShellContainsVertex(const Vertex3 &vert)
{
	if (vert.belongsToShellFacet)
	{
		for (auto &s_facet : shellFacets)
			if (s_facet == vert.belongsToShellFacet)
				return true;
	}
	else if (vert.belongsToShellEdge)
	{
		for (auto &s_edge : shellEdges)
			if (s_edge == vert.belongsToShellEdge)
				return true;
	}
	else if (vert.belongsToShellVertex)
	{
		for (auto &s_edge : shellEdges)
			if (s_edge->vertexes[0] == vert.belongsToShellVertex ||
				s_edge->vertexes[1] == vert.belongsToShellVertex)
				return true;
	}

	return false;
}

void Crystallite3::SetStartFront(
	//const vector<unique_ptr<Vertex3>*> &verts, 
	const vector<unique_ptr<Edge3>*> &edges, 
	const vector<unique_ptr<Facet3>*> &facets)
{
	for (auto &edge : edges)
		if (*edge &&
			ShellContainsVertex(**(*edge)->vertexes[0]) &&
			ShellContainsVertex(**(*edge)->vertexes[1]))
		{
			frontEdges.push_back(
				(new FrontEdge3(edge->get()))
					->GetPtrToUniquePtr());
		}

	for (auto &facet : facets)
		if (*facet &&
			ShellContainsVertex(**(*(*facet)->edges[0])->vertexes[0]) &&
			ShellContainsVertex(**(*(*facet)->edges[0])->vertexes[1]) &&
			ShellContainsVertex(**(*facet)->FindVertexNotIncludedInEdge(**(*facet)->edges[0])))
		{
			frontFacets.push_back(
				(new FrontFacet3(facet->get()))
					->GetPtrToUniquePtr());
		}
}

unique_ptr<FrontFacet3>* Crystallite3::FindFrontFacet(unique_ptr<Facet3>* &facet)
{
	for (auto &f_facet : frontFacets)
		if (*f_facet &&
			(*f_facet)->facet == facet->get())
			return f_facet;

	return nullptr;
}

unique_ptr<FrontEdge3>* Crystallite3::FindFrontEdge(const unique_ptr<Vertex3>* &v0, const unique_ptr<Vertex3>* &v1)
{
	for (auto f_edge : frontEdges)
		if (*f_edge &&
			(((*f_edge)->edge->vertexes[0] == v0 &&
			  (*f_edge)->edge->vertexes[1] == v1) ||
			 ((*f_edge)->edge->vertexes[1] == v0 &&
			  (*f_edge)->edge->vertexes[0] == v1)))
			return f_edge;

	return nullptr;
}

unique_ptr<FrontEdge3>* Crystallite3::FindFrontEdge(unique_ptr<Edge3>* &edge)
{
	for (auto &f_edge : frontEdges)
		if (*f_edge &&
			(*f_edge)->edge == edge->get())
			return f_edge;

	return nullptr;
}

void Crystallite3::AddFrontFacet(unique_ptr<Edge3>* &edge0, unique_ptr<Edge3>* &edge1, unique_ptr<Edge3>* &edge2)
{
	FrontFacet3* new_f_facet = 
		new FrontFacet3(
			new Facet3(
				**edge0, 
				**edge1, 
				**edge2));

	frontFacets.push_back(new_f_facet->GetPtrToUniquePtr());
	innerFacets.push_back(new_f_facet->facet->GetPtrToUniquePtr());
}

const bool Crystallite3::LineSegmentIntersectionCheck(const Vector3 &v0, const Vector3 &v1)
{
	for (auto &f_facet : frontFacets)
		if (*f_facet &&
			Vector3::LineSegmentIntersectTriangle(
				v0,
				v1,
				(*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition(),
				(*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition(),
				(*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition()))
			return true;

	return false;
}

const bool Crystallite3::EdgeIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge)
{
	unique_ptr<Vertex3>* opp_verts[2];
	(*frontEdge)->FindOppositeVertexes(
		frontFacets,
		frontEdges,
		opp_verts[0],
		opp_verts[1]);
	Vector3 opp_verts_poses[2];
	opp_verts_poses[0] = (*opp_verts[0])->GetPosition();
	opp_verts_poses[1] = (*opp_verts[1])->GetPosition();
	Vector3 delta_vec_0th_to_1st = 1e-1 * (opp_verts_poses[1] - opp_verts_poses[0]);

	return LineSegmentIntersectionCheck(
		opp_verts_poses[0] + delta_vec_0th_to_1st,
		opp_verts_poses[1] - delta_vec_0th_to_1st);
}

const bool Crystallite3::FrontSplitCheck(const unique_ptr<FrontEdge3>* frontEdge)
{
	unique_ptr<FrontEdge3>* opp_f_edge = (*frontEdge)->FindOppositeFrontEdge(frontFacets, frontEdges);
	if (!opp_f_edge)
		return false;

	unique_ptr<Vertex3>* opp_f_edge_opp_verts[2];
	(*opp_f_edge)->FindOppositeVertexes(
		frontFacets, 
		frontEdges, 
		opp_f_edge_opp_verts[0], 
		opp_f_edge_opp_verts[1]);

	if (opp_f_edge_opp_verts[0] == (*frontEdge)->edge->vertexes[0] ||
		opp_f_edge_opp_verts[0] == (*frontEdge)->edge->vertexes[1] ||
		opp_f_edge_opp_verts[1] == (*frontEdge)->edge->vertexes[0] ||
		opp_f_edge_opp_verts[1] == (*frontEdge)->edge->vertexes[1])
		return false;

	return true;
}

const bool Crystallite3::FrontContainsOfOnly1Simplex3()
{
	size_t f_facets_num = 0ull;
	for (auto &f_facet : frontFacets)
	{
		if (*f_facet)
		{
			f_facets_num++;
			if (f_facets_num > 4ull)
				return false;
		}
	}

	return true;
}

void Crystallite3::CreateSimplex3WithoutNewVertex(unique_ptr<FrontEdge3>* frontEdge)
{
	unique_ptr<FrontEdge3>* opp_f_edge = (*frontEdge)->FindOppositeFrontEdge(frontFacets, frontEdges);
	if (opp_f_edge)
	{
		unique_ptr<FrontFacet3>* opp_f_facets[2];
		(*opp_f_edge)->FindFrontFacetsAround(
			frontFacets,
			opp_f_facets[0],
			opp_f_facets[1]);

		unique_ptr<FrontFacet3>* main_f_facets[3];
		unique_ptr<Vertex3>* main_vert;
		(*frontEdge)->FindFrontFacetsAround(
			frontFacets,
			main_f_facets[0],
			main_f_facets[1]);
		if ((*opp_f_facets[0])->facet->IsContaining(**(*frontEdge)->edge->vertexes[0]))
		{
			main_f_facets[2] = opp_f_facets[0];
			main_vert = (*frontEdge)->edge->vertexes[0];
		}
		else if ((*opp_f_facets[0])->facet->IsContaining(**(*frontEdge)->edge->vertexes[1]))
		{
			main_f_facets[2] = opp_f_facets[0];
			main_vert = (*frontEdge)->edge->vertexes[1];
		}
		else if ((*opp_f_facets[1])->facet->IsContaining(**(*frontEdge)->edge->vertexes[0]))
		{
			main_f_facets[2] = opp_f_facets[1];
			main_vert = (*frontEdge)->edge->vertexes[0];
		}
		else
		{
			main_f_facets[2] = opp_f_facets[1];
			main_vert = (*frontEdge)->edge->vertexes[1];
		}

		unique_ptr<FrontEdge3>* new_f_facet_f_edges[3];
		unique_ptr<Edge3>* new_facet_edges[3];
		for (int i = 0; i < 3; i++)
		{
			new_facet_edges[i] = (*main_f_facets[i])->facet->FindEdgeNotContainingVertex(**main_vert);
			new_f_facet_f_edges[i] = FindFrontEdge(new_facet_edges[i]);
		}

		//delete frontEdge->release();

		unique_ptr<FrontEdge3>* f_edge_ = nullptr;
		for (auto &f_facet : main_f_facets)
			for (auto &edge : (*f_facet)->facet->edges)
				if ((*edge)->IsContaining(**main_vert) &&
					(f_edge_ = FindFrontEdge(edge)))
					delete f_edge_->release();

		for (auto &f_facet : main_f_facets)
			delete f_facet->release();

		AddFrontFacet(
			new_facet_edges[0],
			new_facet_edges[1],
			new_facet_edges[2]);

		// Maybe add new simplex3...
	}
	else
	{

	}
}

void Crystallite3::ProcessSmallAngles()
{
	if (FrontContainsOfOnly1Simplex3())
	{
		// Do something...

		for (auto &f_facet : frontFacets)
			if (*f_facet)
				delete f_facet->release();
		for (auto &f_edge : frontEdges)
			if (*f_edge)
				delete f_edge->release();
		return;
	}

	for (auto it = frontEdges.begin(); it != frontEdges.end(); it++)
	{
		if (**it)
		{
			if (!FrontSplitCheck(*it) &&
				!EdgeIntersectionCheck(*it) &&
				(**it)->Angle(frontFacets) < DEG_80_IN_RADIANS)
			{
				CreateSimplex3WithoutNewVertex(*it);
				
				if (FrontContainsOfOnly1Simplex3())
				{
					// Do something...

					for (auto &f_facet : frontFacets)
						if (*f_facet)
							delete f_facet->release();
					for (auto &f_edge : frontEdges)
						if (*f_edge)
							delete f_edge->release();
					return;
				}
			}
		}
	}
}

void Crystallite3::TriangulateVolume(const double preferredLength)
{
	ProcessSmallAngles();
	//ProcessMediumAngles();
	//ProcessLargeAngles();
}

Crystallite3::Crystallite3() {}

Crystallite3::~Crystallite3()
{
	//for (auto &simp : innerSimps)
	//{
	//	if (*simp)
	//		delete simp->release();
	//	delete simp;
	//}
	for (auto &facet : innerFacets)
	{
		if (*facet)
			delete facet->release();
		delete facet;
	}
	for (auto &edge : innerEdges)
	{
		if (*edge)
			delete edge->release();
		delete edge;
	}
	for (auto &vert : innerVerts)
	{
		if (*vert)
			delete vert->release();
		delete vert;
	}

	for (auto &ptr : frontFacets)
	{
		if (*ptr)
			delete ptr->release();
		delete ptr;
	}
	for (auto &ptr : frontEdges)
	{
		if (*ptr)
			delete ptr->release();
		delete ptr;
	}
}