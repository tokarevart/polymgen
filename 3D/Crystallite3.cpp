#include "Crystallite3.h"
#include <algorithm>

#define DET(a, b, c, d) \
		(a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))


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
		if (ShellContainsVertex(**(*edge)->vertexes[0]) &&
			ShellContainsVertex(**(*edge)->vertexes[1]))
			frontEdges.push_back(
				(new FrontEdge3(**edge))
					->GetPtrToUniquePtr());

	for (auto &facet : facets)
		if (ShellContainsVertex(**(*(*facet)->edges[0])->vertexes[0]) &&
			ShellContainsVertex(**(*(*facet)->edges[0])->vertexes[1]) &&
			ShellContainsVertex(**(*facet)->FindVertexNotIncludedInEdge(**(*facet)->edges[0])))
			frontFacets.push_back(
				(new FrontFacet3(**facet))
					->GetPtrToUniquePtr());
}

const bool Crystallite3::EdgeIntersectionCheck(const Vector3 &v0, const Vector3 &v1)
{
	for (auto &f_facet : frontFacets)
		if (Vector3::LineSegmentIntersectTriangle(
				v0,
				v1,
				(*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition(),
				(*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition(),
				(*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition()))
			return true;

	return false;
}

void Crystallite3::TriangulateVolume(const double preferredLength)
{
	//ProcessSmallAngles();
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