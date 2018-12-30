#include "Facet3.h"


const Vec3 Facet3::computeCenter()
{
	return 0.3333333333333333 * (
		(*(*edges[0])->vertexes[0])->getPosition() +
		(*(*edges[0])->vertexes[1])->getPosition() +
		(*findVertexNotIncludedInEdge(**edges[0]))->getPosition());
}

const double Facet3::computeQuality()
{
	return (*findShortestEdge())->sqrMagnitude() / (*findLongestEdge())->sqrMagnitude();
}

unique_ptr<Edge3>* Facet3::intersectAlongAnEdge(const Facet3& facet0, const Facet3& facet1)
{
	int inters = 0;
	unique_ptr<Edge3>* res = nullptr;
	for (auto &edge0 : facet0.edges)
		for (auto &edge1 : facet1.edges)
			if (edge0 == edge1)
			{
				inters++;
				res = edge0;
				break;
			}

	return inters == 1 ? res : nullptr;
}

const bool Facet3::intersectsBy(const Point3& origin, const Vec3& dir)
{
	return spatialalgs::isRayIntersectTriangle(
		origin, dir,
		(*(*edges[0])->vertexes[0])->getPosition(), 
		(*(*edges[0])->vertexes[1])->getPosition(), 
		(*findVertexNotIncludedInEdge(**edges[0]))->getPosition());
}

unique_ptr<Vertex3>* Facet3::findVertexNotIncludedInEdge(const Edge3& edge) const
{
	for (auto &facet_edge : edges)
	{
		if (&edge != facet_edge->get())
		{
			if (!edge.contains(**(*facet_edge)->vertexes[0]))
				return (*facet_edge)->vertexes[0];
			else
				return (*facet_edge)->vertexes[1];
		}
	}

	return nullptr;
}

unique_ptr<Edge3>* Facet3::findEdgeNotContainingVertex(const Vertex3& vert) const
{
	for (auto &edge : edges)
		if (!(*edge)->contains(vert))
			return edge;

	return nullptr;
}

unique_ptr<Edge3>* Facet3::findEdge(const Vertex3& vert0, const Vertex3& vert1)
{
	for (auto &edge : edges)
		if (((*edge)->vertexes[0]->get() == &vert0 &&
			 (*edge)->vertexes[1]->get() == &vert1) ||
			((*edge)->vertexes[0]->get() == &vert1 &&
			 (*edge)->vertexes[1]->get() == &vert0))
			return edge;
	
	return nullptr;
}

unique_ptr<Edge3>* Facet3::findShortestEdge()
{
	if ((*edges[0])->sqrMagnitude() < (*edges[1])->sqrMagnitude())
	{
		if ((*edges[0])->sqrMagnitude() < (*edges[2])->sqrMagnitude())
			return edges[0];
		else
			return edges[2];
	}
	else
	{
		if ((*edges[1])->sqrMagnitude() < (*edges[2])->sqrMagnitude())
			return edges[1];
		else
			return edges[2];
	}

	return nullptr;
}

unique_ptr<Edge3>* Facet3::findLongestEdge()
{
	if ((*edges[0])->sqrMagnitude() > (*edges[1])->sqrMagnitude())
	{
		if ((*edges[0])->sqrMagnitude() > (*edges[2])->sqrMagnitude())
			return edges[0];
		else
			return edges[2];
	}
	else
	{
		if ((*edges[1])->sqrMagnitude() > (*edges[2])->sqrMagnitude())
			return edges[1];
		else
			return edges[2];
	}

	return nullptr;
}

const bool Facet3::contains(const Edge3& edge) const
{
	for (auto &edge_ : edges)
		if (edge_->get() == &edge)
			return true;

	return false;
}

const bool Facet3::contains(const Vertex3& vert) const
{
	for (auto &edge : edges)
		if ((*edge)->contains(vert))
			return true;

	return false;
}

Facet3::Facet3() : unique_ptr_helper<Facet3>(this) 
{
	for (auto &edge : edges)
		edge = nullptr;
}

Facet3::Facet3(Edge3& edge0, Edge3& edge1, Edge3& edge2) : unique_ptr_helper<Facet3>(this)
{
	edges[0] = edge0.getPtrToUPtr();
	edges[1] = edge1.getPtrToUPtr();
	edges[2] = edge2.getPtrToUPtr();
}

Facet3::Facet3(Vertex3& vert0, Vertex3& vert1, Vertex3& vert2) : unique_ptr_helper<Facet3>(this)
{
	edges[0] = (new Edge3(vert0, vert1))->getPtrToUPtr();
	edges[1] = (new Edge3(vert1, vert2))->getPtrToUPtr();
	edges[2] = (new Edge3(vert2, vert0))->getPtrToUPtr();
}

Facet3::~Facet3() {}

const Vec3 FrontFacet3::getNormal()
{
	return _normal;
}

const void FrontFacet3::setNormal(const Vec3& vec)
{
	_normal = vec;
}

const Vec3 FrontFacet3::computeNormal(const vector<unique_ptr<FrontFacet3>*>& frontFacets)
{
	Vec3 center = computeCenter();
	Vec3 third_pos = (*facet->findVertexNotIncludedInEdge(**facet->edges[0]))->getPosition();
	Vec3 normal = Vec3::crossProduct(
		(*(*facet->edges[0])->vertexes[0])->getPosition() - third_pos,
		(*(*facet->edges[0])->vertexes[1])->getPosition() - third_pos).normalize();

	// I don't know why it led to better(correct) result.
	Vec3 test_normal_correct_intersect = normal + Vec3(2.1632737147, 1.488313178, -0.71123534278) * 1e-3;

	auto this_uptr = getPtrToUPtr();
	
	int intersects_num = 0;
	for (auto &f_facet : frontFacets)
	{
		if (!*f_facet ||
			f_facet == this_uptr)
			continue;

		if ((*f_facet)->facet->intersectsBy(center, /*normal*/test_normal_correct_intersect))
			intersects_num++;
	}

	return _normal = intersects_num % 2 == 1 ? normal : -normal;
}

const Vec3 FrontFacet3::computeCenter()
{
	return facet->computeCenter();
}

const double FrontFacet3::computeQuality()
{
	return facet->computeQuality();
}

FrontFacet3::FrontFacet3(Facet3* facet) : facet(facet), unique_ptr_helper<FrontFacet3>(this) {}

FrontFacet3::~FrontFacet3() {}