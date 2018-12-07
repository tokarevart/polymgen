#include "Crystallite3.h"
#include <algorithm>
#include <iostream>

#define DET(a, b, c, d) \
		(a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))

#define PI                 3.141592653589793
#define PI_DIV_2           1.5707963267948966
#define DEG_70_IN_RADIANS  1.2217304763960307
#define DEG_80_IN_RADIANS  1.3962634015954636
#define DEG_90_IN_RADIANS  PI_DIV_2
#define DEG_100_IN_RADIANS 1.7453292519943295
#define DEG_150_IN_RADIANS 2.6179938779914943
#define DEG_160_IN_RADIANS 2.7925268031909273
#define DEG_180_IN_RADIANS PI
#define COS_DEG_60   0.5
#define COS_DEG_70   0.342020143325668733
#define COS_DEG_80   0.173648177666930348
#define COS_DEG_90   0.0
#define COS_DEG_100 -0.173648177666930348
#define COS_DEG_120 -0.5
#define COS_DEG_150 -0.866025403784438646
#define COS_DEG_160 -0.939692620785908384
#define COS_DEG_180 -1.0

template <class T>
void Crystallite3::ErasePtrsToNullptr(vector<unique_ptr<T>*> &vec)
{
	size_t real_objs_num = std::count_if(vec.begin(), vec.end(), [](unique_ptr<T>*& ptr) { return *ptr ? true : false; });
	vector<unique_ptr<T>*> buf_vec(real_objs_num);

	//size_t firts_thr_nums = real_objs_num / 2ull;
	//size_t second_thr_nums = real_objs_num - firts_thr_nums;
	//#pragma omp parallel num_threads(2)
	//{
	//#pragma omp single
	//{
	size_t index = 0ull;
	for (size_t i = 0ull; index < real_objs_num; i++)
	{
		if (*vec[i])
		{
			buf_vec[index] = vec[i];
			index++;
		}
		else
		{
			delete vec[i];
		}
	}
	//}
	//#pragma omp single
	//{
	//size_t index2 = 0ull;
	//for (size_t i = vec.size() - 1ull; index2 < second_thr_nums; i--)
	//{
	//	if (*vec[i])
	//	{
	//		buf_vec[real_objs_num - 1ull - index2] = vec[i];
	//		index2++;
	//	}
	//}
	//}
	//}

	vec = std::move(buf_vec);
}

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
	const vector<unique_ptr<Edge3>*> &edges, 
	const vector<unique_ptr<Facet3>*> &facets)
{
	for (auto &edge : edges)
		if (*edge &&
			ShellContainsVertex(**(*edge)->vertexes[0]) &&
			ShellContainsVertex(**(*edge)->vertexes[1]))
		{
			frontEdges.push_back((new FrontEdge3(edge->get()))->GetPtrToUniquePtr());
		}

	for (auto &facet : facets)
		if (*facet &&
			ShellContainsVertex(**(*(*facet)->edges[0])->vertexes[0]) &&
			ShellContainsVertex(**(*(*facet)->edges[0])->vertexes[1]) &&
			ShellContainsVertex(**(*facet)->FindVertexNotIncludedInEdge(**(*facet)->edges[0])))
		{
			frontFacets.push_back((new FrontFacet3(facet->get()))->GetPtrToUniquePtr());
		}
}

ShellEdge3* Crystallite3::FindShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const
{
	for (auto &s_edge : shellEdges)
		if ((s_edge->vertexes[0] == v0  &&
			 s_edge->vertexes[1] == v1) ||
			(s_edge->vertexes[1] == v0  &&
			 s_edge->vertexes[0] == v1))
		{
			return s_edge;
		}

	return nullptr;
}

unique_ptr<FrontFacet3>* Crystallite3::FindFrontFacet(unique_ptr<Facet3>* facet)
{
	for (auto &f_facet : frontFacets)
	{
		if (*f_facet &&
			(*f_facet)->facet == facet->get())
			return f_facet;
	}

	return nullptr;
}

unique_ptr<FrontEdge3>* Crystallite3::FindFrontEdge(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1)
{
	for (auto f_edge : frontEdges)
	{
		if (*f_edge &&
			(((*f_edge)->edge->vertexes[0] == v0 &&
			  (*f_edge)->edge->vertexes[1] == v1) ||
			 ((*f_edge)->edge->vertexes[1] == v0 &&
			  (*f_edge)->edge->vertexes[0] == v1)))
			return f_edge;
	}

	return nullptr;
}

unique_ptr<FrontEdge3>* Crystallite3::FindFrontEdge(unique_ptr<Edge3>* edge)
{
	for (auto &f_edge : frontEdges)
	{
		if (*f_edge &&
			(*f_edge)->edge == edge->get())
			return f_edge;
	}

	return nullptr;
}

unique_ptr<FrontFacet3>* Crystallite3::AddFrontFacet3(unique_ptr<Edge3>* &edge0, unique_ptr<Edge3>* &edge1, unique_ptr<Edge3>* &edge2)
{
	FrontFacet3* new_f_facet = new FrontFacet3(new Facet3(**edge0, **edge1, **edge2));

	frontFacets.push_back(new_f_facet->GetPtrToUniquePtr());
	innerFacets.push_back(new_f_facet->facet->GetPtrToUniquePtr());

	return new_f_facet->GetPtrToUniquePtr();
}

unique_ptr<FrontEdge3>* Crystallite3::AddFrontEdge3(unique_ptr<Vertex3>* &vert0, unique_ptr<Vertex3>* &vert1)
{
	FrontEdge3* new_f_edge = new FrontEdge3(new Edge3(**vert0, **vert1));

	frontEdges.push_back(new_f_edge->GetPtrToUniquePtr());
	innerEdges.push_back(new_f_edge->edge->GetPtrToUniquePtr());

	return new_f_edge->GetPtrToUniquePtr();
}

void Crystallite3::AddFrontEdge3(unique_ptr<FrontEdge3>* &frontEdge)
{
	frontEdges.push_back((*frontEdge)->GetPtrToUniquePtr());
	innerEdges.push_back((*frontEdge)->edge->GetPtrToUniquePtr());
}

const bool Crystallite3::VertexInsideFrontCheck(const Vector3 &v)
{
	int inters_num = 0;
	for (auto &f_facet : frontFacets)
	{
		if (!*f_facet)
			continue;
		
		Vector3 dir = 
			0.3333333333333333 
			* ((*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition()
				+ (*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition()
				+ (*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition()
				- 3.0 * v);

		for (auto &f_facetj : frontFacets)
		{
			if (*f_facetj && 
				Vector3::RayIntersectTriangle(
					v,
					dir,
					(*(*(*f_facetj)->facet->edges[0])->vertexes[0])->GetPosition(),
					(*(*(*f_facetj)->facet->edges[0])->vertexes[1])->GetPosition(),
					(*(*f_facetj)->facet->FindVertexNotIncludedInEdge(**(*f_facetj)->facet->edges[0]))->GetPosition()))
			{
				inters_num++;
			}
		}

		break;
	}

	return inters_num % 2 == 1;
}

const bool Crystallite3::LineSegmentGlobalIntersectionCheck(const Vector3 &v0, const Vector3 &v1)
{
	for (auto &facet : innerFacets)
	{
		if (*facet &&
			Vector3::SegmentIntersectTriangle(
				v0,
				v1,
				(*(*(*facet)->edges[0])->vertexes[0])->GetPosition(),
				(*(*(*facet)->edges[0])->vertexes[1])->GetPosition(),
				(*(*facet)->FindVertexNotIncludedInEdge(**(*facet)->edges[0]))->GetPosition()))
			return true;
	}

	return false;
}

const bool Crystallite3::LineSegmentFrontIntersectionCheck(const Vector3 &v0, const Vector3 &v1)
{
	for (auto &f_facet : frontFacets)
	{
		if (*f_facet &&
			Vector3::SegmentIntersectTriangle(
				v0,
				v1,
				(*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition(),
				(*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition(),
				(*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition()))
			return true;
	}

	return false;
}

const bool Crystallite3::EdgeGlobalIntersectionCheck(const unique_ptr<Edge3>* edge)
{
	Vector3 delta = 1e-3 * (**(*edge)->vertexes[1] - **(*edge)->vertexes[0]);
	Vector3 segment[2];
	segment[0] = (*(*edge)->vertexes[0])->GetPosition() + delta;
	segment[1] = (*(*edge)->vertexes[1])->GetPosition() - delta;

	for (auto &facet : innerFacets)
	{
		if (*facet &&
			!(*facet)->IsContaining(**edge) &&
			Vector3::SegmentIntersectTriangle(
				segment[0],
				segment[1],
				(*(*(*facet)->edges[0])->vertexes[0])->GetPosition(),
				(*(*(*facet)->edges[0])->vertexes[1])->GetPosition(),
				(*(*facet)->FindVertexNotIncludedInEdge(**(*facet)->edges[0]))->GetPosition()))
			return true;
	}

	return false;
}

const bool XOR(const bool &b0, const bool &b1)
{
	return (b0 || b1) && !(b0 && b1);
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
	Vector3 opp_vert0_to_1 = opp_verts_poses[1] - opp_verts_poses[0];
	Vector3 delta_vec_0th_to_1st = 1e-3 * opp_vert0_to_1;

	if (LineSegmentFrontIntersectionCheck(
		opp_verts_poses[0] + delta_vec_0th_to_1st,
		opp_verts_poses[1] - delta_vec_0th_to_1st))
		return true;

	for (auto &f_edge : frontEdges)
	{
		if (*f_edge &&
			(*f_edge)->edge->IsContaining(**opp_verts[0]) &&
			(*f_edge)->edge->IsContaining(**opp_verts[1]))
			return false;
	}

	for (auto &f_edge : frontEdges)
	{
		if (!*f_edge)
			continue;

		bool contains[2];
		contains[0] = (*f_edge)->edge->IsContaining(**opp_verts[0]);
		contains[1] = (*f_edge)->edge->IsContaining(**opp_verts[1]);

		unique_ptr<Vertex3>* vert_buf;
		if (contains[0])
		{
			if ((*f_edge)->edge->vertexes[0] == opp_verts[0])
				vert_buf = (*f_edge)->edge->vertexes[1];
			else
				vert_buf = (*f_edge)->edge->vertexes[0];

			if ((*vert_buf)->GetPosition().DistanceToSegment(
					opp_verts_poses[0],
					opp_verts_poses[1]) < 4e-3 * _preferredLength)
				return true;
		}
		else if (contains[1])
		{
			if ((*f_edge)->edge->vertexes[0] == opp_verts[1])
				vert_buf = (*f_edge)->edge->vertexes[1];
			else
				vert_buf = (*f_edge)->edge->vertexes[0];

			if ((*vert_buf)->GetPosition().DistanceToSegment(
					opp_verts_poses[0],
					opp_verts_poses[1]) < 4e-3 * _preferredLength)
				return true;
		}
		else
		{
			if (Vector3::SegmentsDistance(
					opp_verts_poses[0],
					opp_verts_poses[1],
					(*(*f_edge)->edge->vertexes[0])->GetPosition(),
					(*(*f_edge)->edge->vertexes[1])->GetPosition()) < 4e-3 * _preferredLength)
				return true;
		}
	}

	return false;
}

const bool Crystallite3::FacetIntersectionCheck(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1, const Vector3 &v2)
{
	Vector3 verts_poses[3];
	verts_poses[0] = (*v0)->GetPosition();
	verts_poses[1] = (*v1)->GetPosition();

	for (auto &f_edge : frontEdges)
	{
		if (!*f_edge)
			continue;

		bool contains[3];
		contains[0] = (*f_edge)->edge->IsContaining(**v0);
		contains[1] = (*f_edge)->edge->IsContaining(**v1);

		if (contains[0] && contains[1])
			continue;

		Vector3 first_to_second_delta = 1e-3 * ((*(*f_edge)->edge->vertexes[1])->GetPosition() - (*(*f_edge)->edge->vertexes[0])->GetPosition());
		Vector3 f_edge_verts_poses[2];
		if (contains[0])
		{
			if ((*f_edge)->edge->vertexes[0] == v0)
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition() + first_to_second_delta;
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition();
			}
			else
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition();
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition() - first_to_second_delta;
			}
		}
		else if (contains[1])
		{
			if ((*f_edge)->edge->vertexes[0] == v1)
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition() + first_to_second_delta;
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition();
			}
			else
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition();
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition() - first_to_second_delta;
			}
		}
		else
		{
			f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition();
			f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition();
		}

		if (Vector3::SegmentIntersectTriangle(
			f_edge_verts_poses[0],
			f_edge_verts_poses[1],
			verts_poses[0],
			verts_poses[1],
			v2))
			return true;
	}

	return false;
}

const bool Crystallite3::FacetIntersectionCheck(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1, const unique_ptr<Vertex3>* v2)
{
	Vector3 verts_poses[3];
	verts_poses[0] = (*v0)->GetPosition();
	verts_poses[1] = (*v1)->GetPosition();
	verts_poses[2] = (*v2)->GetPosition();

	for (auto &f_edge : frontEdges)
	{
		if (!*f_edge)
			continue;

		bool contains[3];
		contains[0] = (*f_edge)->edge->IsContaining(**v0);
		contains[1] = (*f_edge)->edge->IsContaining(**v1);
		contains[2] = (*f_edge)->edge->IsContaining(**v2);

		if ((contains[0] && contains[1]) ||
			(contains[0] && contains[2]) ||
			(contains[1] && contains[2]))
			continue;

		Vector3 first_to_second_delta = 1e-3 * ((*(*f_edge)->edge->vertexes[1])->GetPosition() - (*(*f_edge)->edge->vertexes[0])->GetPosition());
		Vector3 f_edge_verts_poses[2];
		if (contains[0])
		{
			if ((*f_edge)->edge->vertexes[0] == v0)
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition() + first_to_second_delta;
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition();
			}
			else
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition();
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition() - first_to_second_delta;
			}
		}
		else if (contains[1])
		{
			if ((*f_edge)->edge->vertexes[0] == v1)
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition() + first_to_second_delta;
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition();
			}
			else
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition();
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition() - first_to_second_delta;
			}
		}
		else if (contains[2])
		{
			if ((*f_edge)->edge->vertexes[0] == v2)
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition() + first_to_second_delta;
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition();
			}
			else
			{
				f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition();
				f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition() - first_to_second_delta;
			}
		}
		else
		{
			f_edge_verts_poses[0] = (*(*f_edge)->edge->vertexes[0])->GetPosition();
			f_edge_verts_poses[1] = (*(*f_edge)->edge->vertexes[1])->GetPosition();
		}

		if (Vector3::SegmentIntersectTriangle(
				f_edge_verts_poses[0],
				f_edge_verts_poses[1],
				verts_poses[0],
				verts_poses[1],
				verts_poses[2]))
			return true;
	}

	return false;
}

const bool Crystallite3::FacetsIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge)
{
	unique_ptr<Vertex3>* opp_verts[2];
	(*frontEdge)->FindOppositeVertexes(frontFacets, frontEdges, opp_verts[0], opp_verts[1]);

	Vector3 opp_verts_poses[2];
	opp_verts_poses[0] = (*opp_verts[0])->GetPosition();
	opp_verts_poses[1] = (*opp_verts[1])->GetPosition();

	//for (auto f_edge : frontEdges)
	//{
	//	if (*f_edge &&
	//		(((*f_edge)->edge->vertexes[0] == opp_verts[0] &&
	//		  (*f_edge)->edge->vertexes[1] == opp_verts[1]) ||
	//		 ((*f_edge)->edge->vertexes[1] == opp_verts[0] &&
	//		  (*f_edge)->edge->vertexes[0] == opp_verts[1])))
	//		return false;
	//}

	if (FacetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[0], opp_verts[1]) ||
		FacetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], opp_verts[1]))
		return true;

	return false;
}

const bool Crystallite3::InsideSimplex3Check(const Vector3 &p0, const Vector3 &p1, const Vector3 &p2, const Vector3 &p3, const Vector3 &vert)
{
	Vector3 vert_to_p0 = p0 - vert;
	Vector3 vert_to_p1 = p1 - vert;
	Vector3 vert_to_p2 = p2 - vert;
	Vector3 vert_to_p3 = p3 - vert;

	double abs_mixed_prods[5];
	abs_mixed_prods[0] = abs(Vector3::MixedProduct(vert_to_p0, vert_to_p2, vert_to_p3));
	abs_mixed_prods[1] = abs(Vector3::MixedProduct(vert_to_p0, vert_to_p1, vert_to_p2));
	abs_mixed_prods[2] = abs(Vector3::MixedProduct(vert_to_p0, vert_to_p1, vert_to_p3));
	abs_mixed_prods[3] = abs(Vector3::MixedProduct(vert_to_p1, vert_to_p2, vert_to_p3));
	abs_mixed_prods[4] = abs(Vector3::MixedProduct(p1 - p0, p2 - p0, p3 - p0));

	/*if (abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4] * 1.0000001)
		return true;
	else
		return false;*/
	return abs_mixed_prods[0] + abs_mixed_prods[1] + abs_mixed_prods[2] + abs_mixed_prods[3] < abs_mixed_prods[4] * 1.000001;
}

const bool Crystallite3::SomeVertexInsidePotentialSimplex3Check(const unique_ptr<FrontEdge3>* frontEdge)
{
	unique_ptr<Vertex3>* opp_verts[2];
	(*frontEdge)->FindOppositeVertexes(frontFacets, frontEdges, opp_verts[0], opp_verts[1]);

	Vector3 points[4];
	points[0] = (*opp_verts[0])->GetPosition();
	points[1] = (*opp_verts[1])->GetPosition();
	points[2] = (*(*frontEdge)->edge->vertexes[0])->GetPosition();
	points[3] = (*(*frontEdge)->edge->vertexes[1])->GetPosition();
	for (auto &vert : innerVerts)
	{
		if (*vert &&
			vert != opp_verts[0] &&
			vert != opp_verts[1] &&
			vert != (*frontEdge)->edge->vertexes[0] &&
			vert != (*frontEdge)->edge->vertexes[1] &&
			InsideSimplex3Check(points[0], points[1], points[2], points[3], (*vert)->GetPosition()))
		{
			return true;
		}
	}

	return false;
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

const bool Crystallite3::ParallelFacetsCheck(const unique_ptr<FrontEdge3>* frontEdge)
{
	unique_ptr<FrontFacet3>* f_facets_around[2];
	(*frontEdge)->FindFrontFacetsAround(frontFacets, f_facets_around[0], f_facets_around[1]);

	unique_ptr<Vertex3>* opp_verts[2];
	opp_verts[0] = (*f_facets_around[0])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);
	opp_verts[1] = (*f_facets_around[1])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);

	Vector3 plane0[2];
	plane0[0] = **opp_verts[0] - **(*frontEdge)->edge->vertexes[0];
	plane0[1] = **opp_verts[1] - **(*frontEdge)->edge->vertexes[0];
	Vector3 normal0 = Vector3::CrossProduct(plane0[0], plane0[1]).Normalize();
	Vector3 plane1[2];
	plane1[0] = **opp_verts[0] - **(*frontEdge)->edge->vertexes[1];
	plane1[1] = **opp_verts[1] - **(*frontEdge)->edge->vertexes[1];
	Vector3 normal1 = Vector3::CrossProduct(plane1[0], plane1[1]).Normalize();

	for (auto &f_facet : frontFacets)
	{
		unique_ptr<Edge3>* inters_reses[2];
		if (*f_facet &&
			f_facet != f_facets_around[0] &&
			f_facet != f_facets_around[1])
		{
			inters_reses[0] = Facet3::IntersectAlongAnEdge(*(*f_facet)->facet, *(*f_facets_around[0])->facet);
			inters_reses[1] = Facet3::IntersectAlongAnEdge(*(*f_facet)->facet, *(*f_facets_around[1])->facet);
			if (!XOR(inters_reses[0], inters_reses[1]))
				continue;

			unique_ptr<Vertex3>* f_facet_to_verts[2];
			f_facet_to_verts[0] = (*(*f_facet)->facet->edges[0])->vertexes[0];
			f_facet_to_verts[1] = (*(*f_facet)->facet->edges[0])->vertexes[1];
			unique_ptr<Vertex3>* f_facet_from_vert = (*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]);

			Vector3 f_plane[2];
			f_plane[0] = **f_facet_to_verts[0] - **f_facet_from_vert;
			f_plane[1] = **f_facet_to_verts[1] - **f_facet_from_vert;
			Vector3 f_normal = Vector3::CrossProduct(f_plane[0], f_plane[1]).Normalize();

			if (abs(abs(Vector3::DotProduct(f_normal, normal0)) - 1.0) < 1e-6 ||
				abs(abs(Vector3::DotProduct(f_normal, normal1)) - 1.0) < 1e-6)
			{
				int i = inters_reses[0] ? 0 : 1;

				Vector3 border_verts[2];
				border_verts[0] = (*(*inters_reses[i])->vertexes[0])->GetPosition();
				border_verts[1] = (*(*inters_reses[i])->vertexes[1])->GetPosition();

				Vector3 main_facet_3rd_vert;
				if ((*inters_reses[i])->IsContaining(**opp_verts[0]))
					main_facet_3rd_vert = (*opp_verts[1])->GetPosition();
				else
					main_facet_3rd_vert = (*opp_verts[0])->GetPosition();

				Vector3 curr_facet_3rd_vert = (*(*f_facet)->facet->FindVertexNotIncludedInEdge(**inters_reses[i]))->GetPosition();

				Vector3 main_facet_cross = Vector3::CrossProduct(main_facet_3rd_vert - border_verts[0], main_facet_3rd_vert - border_verts[1]);
				Vector3 curr_facet_cross = Vector3::CrossProduct(curr_facet_3rd_vert - border_verts[0], curr_facet_3rd_vert - border_verts[1]);
				if (Vector3::DotProduct(main_facet_cross, curr_facet_cross) > 0.0)
					return true;
			}
		}
	}

	return false;
}

const bool Crystallite3::FrontContainsOfOnly1FacetOrEmpty()
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

unique_ptr<FrontEdge3>* Crystallite3::CurrentFrontEdge(double maxExCos)
{
	double curr_max_excos = -2.0;
	unique_ptr<FrontEdge3>* curr_max_f_edge = nullptr;
	for (size_t i = 0ull; i < frontEdges.size(); i++)
	{
		if (!*frontEdges[i])
			continue;

		double curr_excos = (*frontEdges[i])->AngleExCos(frontFacets);
		if (curr_excos > curr_max_excos &&
			curr_excos < maxExCos)
		{
			curr_max_excos = curr_excos;
			curr_max_f_edge = frontEdges[i];
		}
	}

	return curr_max_f_edge;
}

const bool Crystallite3::ExhaustWithoutNewVertexPriorityPredicate(unique_ptr<FrontEdge3>* currentFrontEdge)
{
	if ((*currentFrontEdge)->AngleExCos(frontFacets) > COS_DEG_70)
		return true;

	unique_ptr<Vertex3>* opp_verts[2];
	(*currentFrontEdge)->FindOppositeVertexes(frontFacets, frontEdges, opp_verts[0], opp_verts[1]);
	if (FindFrontEdge(opp_verts[0], opp_verts[1]) ||
		((*currentFrontEdge)->AngleExCos(frontFacets) < COS_DEG_70 &&
		 (*currentFrontEdge)->AngleExCos(frontFacets) > COS_DEG_100 &&
		 (**opp_verts[1] - **opp_verts[0]).SqrMagnitude() <= _preferredLength * _preferredLength))
		return true;

	return false;
}

const bool Crystallite3::ExhaustWithNewVertexPriorityPredicate(unique_ptr<FrontEdge3>* currentFrontEdge)
{
	if ((*currentFrontEdge)->AngleExCos(frontFacets) < COS_DEG_120)
		return true;

	return false;
}

void Crystallite3::ExhaustWithoutNewVertexOppositeEdgeExists(unique_ptr<FrontEdge3>* frontEdge, unique_ptr<FrontEdge3>* oppositeEdge)
{
	unique_ptr<FrontEdge3>* opp_f_edge = oppositeEdge;

	unique_ptr<FrontFacet3>* opp_f_facets[2];
	(*opp_f_edge)->FindFrontFacetsAround(
		frontFacets,
		opp_f_facets[0], opp_f_facets[1]);

	unique_ptr<FrontFacet3>* main_f_facets[3];
	unique_ptr<Vertex3>* main_vert;
	(*frontEdge)->FindFrontFacetsAround(
		frontFacets,
		main_f_facets[0], main_f_facets[1]);
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

	unique_ptr<FrontEdge3>* new_f_facet_f_edge;
	unique_ptr<Edge3>* new_facet_edges[3];
	for (int i = 0; i < 3; i++)
	{
		new_facet_edges[i] = (*main_f_facets[i])->facet->FindEdgeNotContainingVertex(**main_vert);
		new_f_facet_f_edge = FindFrontEdge(new_facet_edges[i]);
		(*new_f_facet_f_edge)->needProcessing = true;
	}

	AddFrontFacet3(
		new_facet_edges[0],
		new_facet_edges[1],
		new_facet_edges[2]);

	innerSimps.push_back((new Simplex3(
		**(*frontEdge)->edge->vertexes[0],
		**(*frontEdge)->edge->vertexes[1],
		**(*opp_f_edge)->edge->vertexes[0],
		**(*opp_f_edge)->edge->vertexes[1]))->GetPtrToUniquePtr());

	unique_ptr<FrontEdge3>* f_edge_ = nullptr;
	for (auto &f_facet : main_f_facets)
		for (auto &edge : (*f_facet)->facet->edges)
			if ((*edge)->IsContaining(**main_vert) &&
				(f_edge_ = FindFrontEdge(edge)))
				delete f_edge_->release();

	for (auto &f_facet : main_f_facets)
		delete f_facet->release();
}

void Crystallite3::ExhaustWithoutNewVertexOppositeEdgeDontExists(unique_ptr<FrontEdge3>* frontEdge)
{
	unique_ptr<FrontEdge3>* opp_f_edge = nullptr;

	unique_ptr<FrontFacet3>* f_facets_around[2];
	(*frontEdge)->FindFrontFacetsAround(
		frontFacets,
		f_facets_around[0], f_facets_around[1]);

	unique_ptr<Vertex3>* opp_verts[2];
	opp_verts[0] = (*f_facets_around[0])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);
	opp_verts[1] = (*f_facets_around[1])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);

	opp_f_edge = (new FrontEdge3(new Edge3(**opp_verts[0], **opp_verts[1])))->GetPtrToUniquePtr();
	AddFrontEdge3(opp_f_edge);

	unique_ptr<FrontEdge3>* new_f_facet_f_edge;
	unique_ptr<Edge3>* new_facet_edges[3];
	new_facet_edges[2] = (*opp_f_edge)->edge->GetPtrToUniquePtr();

	for (auto &edge : (*f_facets_around[0])->facet->edges)
		if ((*edge)->IsContaining(**(*frontEdge)->edge->vertexes[0]) &&
			((*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[0]) ||
			(*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[1])))
		{
			new_facet_edges[0] = edge;
			new_f_facet_f_edge = FindFrontEdge(edge);
			if ((*new_f_facet_f_edge)->needProcessing == false)
				(*new_f_facet_f_edge)->needProcessing = true;
			break;
		}
	for (auto &edge : (*f_facets_around[1])->facet->edges)
		if ((*edge)->IsContaining(**(*frontEdge)->edge->vertexes[0]) &&
			((*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[0]) ||
			(*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[1])))
		{
			new_facet_edges[1] = edge;
			new_f_facet_f_edge = FindFrontEdge(edge);
			if ((*new_f_facet_f_edge)->needProcessing == false)
				(*new_f_facet_f_edge)->needProcessing = true;
			break;
		}
	AddFrontFacet3(
		new_facet_edges[0],
		new_facet_edges[1],
		new_facet_edges[2]);

	for (auto &edge : (*f_facets_around[0])->facet->edges)
		if ((*edge)->IsContaining(**(*frontEdge)->edge->vertexes[1]) &&
			((*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[0]) ||
			(*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[1])))
		{
			new_facet_edges[0] = edge;
			new_f_facet_f_edge = FindFrontEdge(edge);
			if ((*new_f_facet_f_edge)->needProcessing == false)
				(*new_f_facet_f_edge)->needProcessing = true;
			break;
		}
	for (auto &edge : (*f_facets_around[1])->facet->edges)
		if ((*edge)->IsContaining(**(*frontEdge)->edge->vertexes[1]) &&
			((*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[0]) ||
			(*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[1])))
		{
			new_facet_edges[1] = edge;
			new_f_facet_f_edge = FindFrontEdge(edge);
			if ((*new_f_facet_f_edge)->needProcessing == false)
				(*new_f_facet_f_edge)->needProcessing = true;
			break;
		}
	AddFrontFacet3(
		new_facet_edges[0],
		new_facet_edges[1],
		new_facet_edges[2]);

	innerSimps.push_back((new Simplex3(
		**(*frontEdge)->edge->vertexes[0],
		**(*frontEdge)->edge->vertexes[1],
		**opp_verts[0],
		**opp_verts[1]))->GetPtrToUniquePtr());

	delete f_facets_around[0]->release();
	delete f_facets_around[1]->release();
	delete frontEdge->release();
}

void Crystallite3::ExhaustWithoutNewVertex(unique_ptr<FrontEdge3>* frontEdge, const bool oppositeEdgeExistence, unique_ptr<FrontEdge3>* oppositeEdge)
{
	unique_ptr<FrontEdge3>* opp_f_edge = nullptr;
	if (oppositeEdgeExistence && oppositeEdge)
		opp_f_edge = oppositeEdge;
	else if (oppositeEdgeExistence)
		opp_f_edge = (*frontEdge)->FindOppositeFrontEdge(frontFacets, frontEdges);

	if ((oppositeEdge && oppositeEdgeExistence) ||
		opp_f_edge)
	{
		ExhaustWithoutNewVertexOppositeEdgeExists(frontEdge, opp_f_edge);
	}
	else
	{
		ExhaustWithoutNewVertexOppositeEdgeDontExists(frontEdge);
	}
}

const bool Crystallite3::NewVertexPosition_OLD(unique_ptr<FrontEdge3>* frontEdge, Vector3& out_pos)
{
	Vector3 f_edge_verts_poses[2];
	f_edge_verts_poses[0] = (*(*frontEdge)->edge->vertexes[0])->GetPosition();
	f_edge_verts_poses[1] = (*(*frontEdge)->edge->vertexes[1])->GetPosition();

	Vector3 orig = 0.5 * (f_edge_verts_poses[0] + f_edge_verts_poses[1]);

	unique_ptr<FrontFacet3>* f_facets[2];
	(*frontEdge)->FindFrontFacetsAround(frontFacets, f_facets[0], f_facets[1]);

	unique_ptr<Vertex3>* opp_verts[2];
	opp_verts[0] = (*f_facets[0])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);
	opp_verts[1] = (*f_facets[1])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);

	Vector3 opp_verts_poses[2];
	opp_verts_poses[0] = (*opp_verts[0])->GetPosition();
	opp_verts_poses[1] = (*opp_verts[1])->GetPosition();

	Vector3 facets_inner_vecs[2];
	facets_inner_vecs[0] = opp_verts_poses[0] - orig;
	facets_inner_vecs[1] = opp_verts_poses[1] - orig;

	Vector3	dir = (facets_inner_vecs[0] + facets_inner_vecs[1]).Normalize() * _preferredLength * 0.866; // sind(60) == 0.8660...
	Vector3 orig_plus_ex_dir = orig + dir * 2.0;
	while (
		LineSegmentFrontIntersectionCheck(f_edge_verts_poses[0] + 1e-3 * (orig_plus_ex_dir - f_edge_verts_poses[0]), orig_plus_ex_dir) ||
		LineSegmentFrontIntersectionCheck(f_edge_verts_poses[1] + 1e-3 * (orig_plus_ex_dir - f_edge_verts_poses[1]), orig_plus_ex_dir) ||
		LineSegmentFrontIntersectionCheck(opp_verts_poses[0] + 1e-3 * (orig_plus_ex_dir - opp_verts_poses[0]), orig_plus_ex_dir) ||
		LineSegmentFrontIntersectionCheck(opp_verts_poses[1] + 1e-3 * (orig_plus_ex_dir - opp_verts_poses[1]), orig_plus_ex_dir) ||
		!VertexInsideFrontCheck(orig_plus_ex_dir) ||
		FacetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[0], orig_plus_ex_dir) ||
		FacetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], orig_plus_ex_dir) ||
		FacetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[1], orig_plus_ex_dir) ||
		FacetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], orig_plus_ex_dir)
		//!VertexInsideFrontCheck(orig + dir) ||
		/*LineSegmentGlobalIntersectionCheck(f_edge_verts_poses[0] + 1e-3 * (orig + dir - f_edge_verts_poses[0]), orig + dir) ||
		LineSegmentGlobalIntersectionCheck(f_edge_verts_poses[1] + 1e-3 * (orig + dir - f_edge_verts_poses[1]), orig + dir) ||
		LineSegmentGlobalIntersectionCheck(f_facets_far_verts_poses[0] + 1e-3 * (orig + dir - f_facets_far_verts_poses[0]), orig + dir) ||
		LineSegmentFrontIntersectionCheck(f_facets_far_verts_poses[1] + 1e-3 * (orig + dir - f_facets_far_verts_poses[1]), orig + dir)*/)
	{
		dir *= 0.75;
		orig_plus_ex_dir = orig + dir * 2.0;

		if ((frontEdges.size() < 10ull && dir.Magnitude() < _preferredLength * 0.01) ||
			dir.Magnitude() < _preferredLength * 0.05)
			return false;
	}

	/*if (dir.Magnitude() < _preferredLength * 0.03)
		return false;*/

	/*if (FacetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[0], orig + dir) ||
		FacetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], orig + dir) ||
		FacetIntersectionCheck((*frontEdge)->edge->vertexes[0], opp_verts[1], orig + dir) ||
		FacetIntersectionCheck((*frontEdge)->edge->vertexes[1], opp_verts[0], orig + dir))
		return false;*/

	/*for (auto &vert : innerVerts)
	{
		if (!*vert ||
			vert == (*frontEdge)->edge->vertexes[0] ||
			vert == (*frontEdge)->edge->vertexes[1] ||
			vert == opp_verts[0] ||
			vert == opp_verts[1])
			continue;

		if (InsideSimplex3Check(
				(*(*frontEdge)->edge->vertexes[0])->GetPosition(),
				(*(*frontEdge)->edge->vertexes[1])->GetPosition(),
				(*opp_verts[0])->GetPosition(),
				orig_plus_ex_dir,
				(*vert)->GetPosition()) ||
			InsideSimplex3Check(
				(*(*frontEdge)->edge->vertexes[0])->GetPosition(),
				(*(*frontEdge)->edge->vertexes[1])->GetPosition(),
				(*opp_verts[1])->GetPosition(),
				orig_plus_ex_dir,
				(*vert)->GetPosition()))
			return false;
	}*/

	out_pos = orig + dir;
	return true;
}

void Crystallite3::ExhaustWithNewVertex_OLD(unique_ptr<FrontEdge3>* frontEdge, Vector3 vertPos)
{
	unique_ptr<FrontFacet3>* f_facets_around[2];
	(*frontEdge)->FindFrontFacetsAround(
		frontFacets,
		f_facets_around[0], f_facets_around[1]);

	unique_ptr<Vertex3>* opp_verts[2];
	opp_verts[0] = (*f_facets_around[0])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);
	opp_verts[1] = (*f_facets_around[1])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);

	unique_ptr<Vertex3>* new_vert = (new Vertex3(vertPos))->GetPtrToUniquePtr();
	innerVerts.push_back(new_vert);

	unique_ptr<Edge3>* new_facet_edges[9];
	new_facet_edges[0] = (*frontEdge)->edge->GetPtrToUniquePtr();
	new_facet_edges[1] = (*AddFrontEdge3((*frontEdge)->edge->vertexes[0], new_vert))->edge->GetPtrToUniquePtr();
	new_facet_edges[2] = (*AddFrontEdge3((*frontEdge)->edge->vertexes[1], new_vert))->edge->GetPtrToUniquePtr();
	new_facet_edges[3] = (*f_facets_around[0])->facet->FindEdge(**(*frontEdge)->edge->vertexes[0], **opp_verts[0]);
	new_facet_edges[4] = (*AddFrontEdge3(opp_verts[0], new_vert))->edge->GetPtrToUniquePtr();
	new_facet_edges[5] = (*f_facets_around[0])->facet->FindEdge(**(*frontEdge)->edge->vertexes[1], **opp_verts[0]);
	new_facet_edges[6] = (*f_facets_around[1])->facet->FindEdge(**(*frontEdge)->edge->vertexes[0], **opp_verts[1]);
	new_facet_edges[7] = (*AddFrontEdge3(opp_verts[1], new_vert))->edge->GetPtrToUniquePtr();
	new_facet_edges[8] = (*f_facets_around[1])->facet->FindEdge(**(*frontEdge)->edge->vertexes[1], **opp_verts[1]);

	(*FindFrontEdge(new_facet_edges[3]))->needProcessing = true;
	(*FindFrontEdge(new_facet_edges[5]))->needProcessing = true;
	(*FindFrontEdge(new_facet_edges[6]))->needProcessing = true;
	(*FindFrontEdge(new_facet_edges[8]))->needProcessing = true;

	innerFacets.push_back((new Facet3(
		**new_facet_edges[0], 
		**new_facet_edges[1], 
		**new_facet_edges[2]))->GetPtrToUniquePtr());
	AddFrontFacet3(
		new_facet_edges[1],
		new_facet_edges[3],
		new_facet_edges[4]);
	AddFrontFacet3(
		new_facet_edges[2],
		new_facet_edges[4],
		new_facet_edges[5]);
	AddFrontFacet3(
		new_facet_edges[1],
		new_facet_edges[6],
		new_facet_edges[7]);
	AddFrontFacet3(
		new_facet_edges[2],
		new_facet_edges[8],
		new_facet_edges[7]);

	innerSimps.push_back((new Simplex3(
		**(*new_facet_edges[0])->vertexes[0],
		**(*new_facet_edges[0])->vertexes[1],
		**(*new_facet_edges[4])->vertexes[0],
		**(*new_facet_edges[4])->vertexes[1]))->GetPtrToUniquePtr());
	innerSimps.push_back((new Simplex3(
		**(*new_facet_edges[0])->vertexes[0],
		**(*new_facet_edges[0])->vertexes[1],
		**(*new_facet_edges[7])->vertexes[0],
		**(*new_facet_edges[7])->vertexes[1]))->GetPtrToUniquePtr());

	delete f_facets_around[0]->release();
	delete f_facets_around[1]->release();
	delete frontEdge->release();
}

template<class T>
size_t PtrsToNullptrNumber(vector<unique_ptr<T>*> &vec)
{
	return std::count_if(vec.begin(), vec.end(), [](unique_ptr<T>* &ptr) { return !*ptr ? true : false; });
}

//const bool Crystallite3::ProcessAngles(double lessThan, double previousAngles, bool addNewVert, bool isFirst, bool isFinal, Polycrystal3* polycr)
//{
//	if (FrontContainsOfOnly1Simplex3OrEmpty())
//	{
//		// Do something...
//
//		for (auto &f_facet : frontFacets)
//		{
//			if (*f_facet)
//				delete f_facet->release();
//		}
//		for (auto &f_edge : frontEdges)
//		{
//			if (*f_edge)
//				delete f_edge->release();
//		}
//		return true;
//	}
//
//	for (size_t i = 0ull; i < frontEdges.size();)
//	{
//		Vector3 new_vert_pos;
//		if (isFirst)
//		{
//			if (!*frontEdges[i] ||
//				(*frontEdges[i])->Angle(frontFacets) >= lessThan ||
//				ParallelFacetsCheck(frontEdges[i]) ||
//				EdgeIntersectionCheck(frontEdges[i]) ||
//				FacetsIntersectionCheck(frontEdges[i]) ||
//				SomeVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
//				FrontSplitCheck(frontEdges[i]))
//			{
//				i++;
//				continue;
//			}
//		}
//		else
//		{
//			if (!*frontEdges[i] ||
//				(*frontEdges[i])->Angle(frontFacets) >= lessThan ||
//				EdgeIntersectionCheck(frontEdges[i]) ||
//				ParallelFacetsCheck(frontEdges[i]) ||
//				!NewVertexPosition_OLD(frontEdges[i], new_vert_pos))
//			{
//				i++;
//				continue;
//			}
//		}
//
//		if (addNewVert)
//			ExhaustWithNewVertex_OLD(frontEdges[i], new_vert_pos);
//		else
//			ExhaustWithoutNewVertex(frontEdges[i]);
//		i = 0ull;
//
//		if (GlobalIntersectionCheck())
//			std::cout << "Bad!";
//		polycr->OutputData();
//
//		if (ProcessAngles(previousAngles  polycr))
//			return true;
//		if (GlobalIntersectionCheck())
//			std::cout << "Bad!";
//		polycr->OutputData();
//
//		if (FrontContainsOfOnly1Simplex3OrEmpty())
//		{
//			// Do something...
//
//			for (auto &f_facet : frontFacets)
//			{
//				if (*f_facet)
//					delete f_facet->release();
//			}
//			for (auto &f_edge : frontEdges)
//			{
//				if (*f_edge)
//					delete f_edge->release();
//			}
//			return true;
//		}
//
//		if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 3ull)
//			ErasePtrsToNullptr(frontEdges);
//		if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 3ull)
//			ErasePtrsToNullptr(frontFacets);
//	}
//
//	return false;
//}

const bool Crystallite3::ProcessVerySmallAngles(Polycrystal3* polycr)
{
	if (FrontContainsOfOnly1FacetOrEmpty())
	{
		// Do something...

		for (auto &f_facet : frontFacets)
		{
			if (*f_facet)
				delete f_facet->release();
		}
		for (auto &f_edge : frontEdges)
		{
			if (*f_edge)
				delete f_edge->release();
		}
		return true;
	}

	for (size_t i = 0ull; i < frontEdges.size();)
	{
		if (!*frontEdges[i])
		{
			i++;
			continue;
		}

		if ((*frontEdges[i])->AngleExCos(frontFacets) <= COS_DEG_80 ||
			ParallelFacetsCheck(frontEdges[i]) ||
			EdgeIntersectionCheck(frontEdges[i]) ||
			FacetsIntersectionCheck(frontEdges[i]) ||
			SomeVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
			FrontSplitCheck(frontEdges[i]))
		{
			i++;
			continue;
		}

		/*if (!*frontEdges[i] || 
			(*frontEdges[i])->Angle(frontFacets) >= DEG_80_IN_RADIANS ||
			ParallelFacetsCheck(frontEdges[i]) ||
			EdgeIntersectionCheck(frontEdges[i]) ||
			FacetsIntersectionCheck(frontEdges[i]) ||
			SomeVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
			FrontSplitCheck(frontEdges[i]))
		{
			i++;
			continue;
		}*/

		ExhaustWithoutNewVertex(frontEdges[i]);
		i = 0ull;

		/*if (GlobalIntersectionCheck())
		std::cout << "Bad!";*/
		//polycr->OutputData();

		if (FrontContainsOfOnly1FacetOrEmpty())
		{
			// Do something...

			for (auto &f_facet : frontFacets)
			{
				if (*f_facet)
					delete f_facet->release();
			}
			for (auto &f_edge : frontEdges)
			{
				if (*f_edge)
					delete f_edge->release();
			}
			return true;
		}

		if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
			ErasePtrsToNullptr(frontEdges);
		if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
			ErasePtrsToNullptr(frontFacets);
	}

	return false;
}

const bool Crystallite3::ProcessSmallAngles(Polycrystal3* polycr)
{
	if (FrontContainsOfOnly1FacetOrEmpty())
	{
		// Do something...

		for (auto &f_facet : frontFacets)
		{
			if (*f_facet)
				delete f_facet->release();
		}
		for (auto &f_edge : frontEdges)
		{
			if (*f_edge)
				delete f_edge->release();
		}
		return true;
	}

	for (size_t i = 0ull; i < frontEdges.size();)
	{
		if (!*frontEdges[i])
		{
			i++;
			continue;
		}

		if ((*frontEdges[i])->AngleExCos(frontFacets) <= COS_DEG_100 ||
			ParallelFacetsCheck(frontEdges[i]) ||
			EdgeIntersectionCheck(frontEdges[i]) ||
			FacetsIntersectionCheck(frontEdges[i]) ||
			SomeVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
			FrontSplitCheck(frontEdges[i]))
		{
			i++;
			continue;
		}

		/*if (!*frontEdges[i] ||
			(*frontEdges[i])->Angle(frontFacets) >= DEG_100_IN_RADIANS ||
			ParallelFacetsCheck(frontEdges[i]) ||
			EdgeIntersectionCheck(frontEdges[i]) ||
			FacetsIntersectionCheck(frontEdges[i]) ||
			SomeVertexInsidePotentialSimplex3Check(frontEdges[i]) ||
			FrontSplitCheck(frontEdges[i]))
		{
			i++;
			continue;
		}*/

		ExhaustWithoutNewVertex(frontEdges[i]);
		i = 0ull;

		/*if (GlobalIntersectionCheck())
			std::cout << "Bad!";*/
		//polycr->OutputData();

		if (FrontContainsOfOnly1FacetOrEmpty())
		{
			// Do something...

			for (auto &f_facet : frontFacets)
			{
				if (*f_facet)
					delete f_facet->release();
			}
			for (auto &f_edge : frontEdges)
			{
				if (*f_edge)
					delete f_edge->release();
			}
			return true;
		}

		if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
			ErasePtrsToNullptr(frontEdges);
		if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
			ErasePtrsToNullptr(frontFacets);
	}

	return false;
}

const bool Crystallite3::ProcessMediumAngles(Polycrystal3* polycr)
{
	if (FrontContainsOfOnly1FacetOrEmpty())
	{
		// Do something...

		for (auto &f_facet : frontFacets)
		{
			if (*f_facet)
				delete f_facet->release();
		}
		for (auto &f_edge : frontEdges)
		{
			if (*f_edge)
				delete f_edge->release();
		}
		return true;
	}

	for (size_t i = 0ull; i < frontEdges.size();)
	{
		if (!*frontEdges[i])
		{
			i++;
			continue;
		}

		Vector3 new_vert_pos;
		if ((*frontEdges[i])->AngleExCos(frontFacets) <= COS_DEG_150 ||
			EdgeIntersectionCheck(frontEdges[i]) ||
			ParallelFacetsCheck(frontEdges[i]) ||
			!NewVertexPosition_OLD(frontEdges[i], new_vert_pos))
		{
			i++;
			continue;
		}

		/*Vector3 new_vert_pos;
		if (!*frontEdges[i] ||
			(*frontEdges[i])->Angle(frontFacets) >= DEG_150_IN_RADIANS ||
			EdgeIntersectionCheck(frontEdges[i]) ||
			ParallelFacetsCheck(frontEdges[i]) ||
			!NewVertexPosition_OLD(frontEdges[i], new_vert_pos))
		{
			i++;
			continue;
		}*/

		ExhaustWithNewVertex_OLD(frontEdges[i], new_vert_pos);
		if (ProcessVerySmallAngles(polycr))
			return true;
		if (ProcessSmallAngles(polycr))
			return true;
		//polycr->OutputData();
		i = 0ull;

		if (FrontContainsOfOnly1FacetOrEmpty())
		{
			// Do something...

			for (auto &f_facet : frontFacets)
			{
				if (*f_facet)
					delete f_facet->release();
			}
			for (auto &f_edge : frontEdges)
			{
				if (*f_edge)
					delete f_edge->release();
			}

			return true;
		}

		if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
			ErasePtrsToNullptr(frontEdges);
		if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
			ErasePtrsToNullptr(frontFacets);
	}

	return false;
}

void Crystallite3::ProcessLargeAngles(Polycrystal3* polycr)
{
	if (FrontContainsOfOnly1FacetOrEmpty())
	{
		// Do something...

		for (auto &f_facet : frontFacets)
		{
			if (*f_facet)
				delete f_facet->release();
		}
		for (auto &f_edge : frontEdges)
		{
			if (*f_edge)
				delete f_edge->release();
		}
		return;
	}
	
	for (size_t i = 0ull; i < frontEdges.size();)
	{
		if (!*frontEdges[i])
		{
			i++;
			continue;
		}

		Vector3 new_vert_pos;
		if ((*frontEdges[i])->AngleExCos(frontFacets) <= COS_DEG_180 ||
			EdgeIntersectionCheck(frontEdges[i]) ||
			ParallelFacetsCheck(frontEdges[i]) ||
			!NewVertexPosition_OLD(frontEdges[i], new_vert_pos))
		{
			i++;
			continue;
		}

		ExhaustWithNewVertex_OLD(frontEdges[i], new_vert_pos);
		/*if (GlobalIntersectionCheck())
			std::cout << "Bad!";*/
		//polycr->OutputData();
		if (ProcessMediumAngles(polycr))
			return;
		/*if (GlobalIntersectionCheck())
			std::cout << "Bad!";*/
		//polycr->OutputData();
		i = 0ull;

		if (FrontContainsOfOnly1FacetOrEmpty())
		{
			// Do something...

			for (auto &f_facet : frontFacets)
			{
				if (*f_facet)
					delete f_facet->release();
			}
			for (auto &f_edge : frontEdges)
			{
				if (*f_edge)
					delete f_edge->release();
			}
			return;
		}

		if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
			ErasePtrsToNullptr(frontEdges);
		if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
			ErasePtrsToNullptr(frontFacets);
	}

	polycr->OutputData();
	if (!FrontContainsOfOnly1FacetOrEmpty())
		throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");
}

void Crystallite3::ProcessAngles(Polycrystal3* polycr)
{
	if (FrontContainsOfOnly1FacetOrEmpty())
	{
		throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");

		for (auto &f_facet : frontFacets)
		{
			if (*f_facet)
				delete f_facet->release();
		}
		for (auto &f_edge : frontEdges)
		{
			if (*f_edge)
				delete f_edge->release();
		}
		return;
	}

	double max_excos = 2.0;
	for (unique_ptr<FrontEdge3>* curr_f_edge = CurrentFrontEdge(max_excos); ; curr_f_edge = CurrentFrontEdge(max_excos))
	{
		if (!curr_f_edge)
		{
			polycr->OutputData();
			throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");
		}

		//double curr_f_edge_magn = (*curr_f_edge)->edge->Magnitude();
		//if (curr_f_edge_magn < 0.2 * _preferredLength ||
		//	curr_f_edge_magn > 2.5 * _preferredLength)
		//{
		//	//for (int i = 0; i < 3; i++)
		//	//{
		//		SmoothAroundFrontVertex((*curr_f_edge)->edge->vertexes[0]);
		//		SmoothAroundFrontVertex((*curr_f_edge)->edge->vertexes[1]);
		//	//}
		//	for (auto f_edge : frontEdges)
		//	{
		//		if (!*f_edge)
		//			continue;
		//
		//		if ((*f_edge)->edge->IsContaining(**(*curr_f_edge)->edge->vertexes[0]) ||
		//			(*f_edge)->edge->IsContaining(**(*curr_f_edge)->edge->vertexes[1]))
		//		{
		//			(*f_edge)->needProcessing = true;
		//		}
		//	}
		//
		//	max_excos = 2.0;
		//	continue;
		//}
		//
		//if ((*curr_f_edge)->AngleExCos(frontFacets) > COS_DEG_100)
		//{
		//	if (ParallelFacetsCheck(curr_f_edge) ||
		//		EdgeIntersectionCheck(curr_f_edge) ||
		//		FacetsIntersectionCheck(curr_f_edge) ||
		//		SomeVertexInsidePotentialSimplex3Check(curr_f_edge) ||
		//		FrontSplitCheck(curr_f_edge))
		//	{
		//		//if ((*curr_f_edge)->AngleExCos(frontFacets) < COS_DEG_60 &&
		//		//	!FrontSplitCheck(curr_f_edge))
		//		//{
		//		//	Vector3 new_vert_pos;
		//		//	if (EdgeIntersectionCheck(curr_f_edge) ||
		//		//		ParallelFacetsCheck(curr_f_edge) ||
		//		//		!NewVertexPosition_OLD(curr_f_edge, new_vert_pos)) // Experemental
		//		//	{
		//		//		max_excos = (*curr_f_edge)->AngleExCos(frontFacets);
		//		//		continue;
		//		//	}
		//		//
		//		//	ExhaustWithNewVertex_OLD(curr_f_edge, new_vert_pos);
		//		//	max_excos = 2.0;
		//		//}
		//		//else 
		//		//{
		//		//	max_excos = (*curr_f_edge)->AngleExCos(frontFacets);
		//		//}
		//		//continue;
		//
		//		max_excos = (*curr_f_edge)->AngleExCos(frontFacets);
		//		continue;
		//	}
		//
		//	ExhaustWithoutNewVertex(curr_f_edge);
		//	max_excos = 2.0;
		//}
		//else
		//{
		//	Vector3 new_vert_pos;
		//	if (EdgeIntersectionCheck(curr_f_edge) ||
		//		ParallelFacetsCheck(curr_f_edge) ||
		//		!NewVertexPosition_OLD(curr_f_edge, new_vert_pos))
		//	{
		//		max_excos = (*curr_f_edge)->AngleExCos(frontFacets);
		//		continue;
		//	}
		//
		//	ExhaustWithNewVertex_OLD(curr_f_edge, new_vert_pos);
		//	max_excos = 2.0;
		//}
		
		polycr->OutputData();

		if (FrontContainsOfOnly1FacetOrEmpty())
		{
			// Do something...

			for (auto &f_facet : frontFacets)
			{
				if (*f_facet)
					delete f_facet->release();
			}
			for (auto &f_edge : frontEdges)
			{
				if (*f_edge)
					delete f_edge->release();
			}
			return;
		}

		if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 2ull)
			ErasePtrsToNullptr(frontEdges);
		if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 2ull)
			ErasePtrsToNullptr(frontFacets);
	}

	polycr->OutputData();
	if (!FrontContainsOfOnly1FacetOrEmpty())
		throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");
}

void Crystallite3::TriangulateVolume(const double preferredLength, Polycrystal3* polycr)
{
	_preferredLength = preferredLength;
	////polycr->OutputData();
	//if (ProcessVerySmallAngles(polycr))
	//{
	//	if (GlobalIntersectionCheck())
	//		std::cout << "Bad!";
	//	return;
	//}
	////polycr->OutputData();
	//if (ProcessSmallAngles(polycr))
	//{
	//	if (GlobalIntersectionCheck())
	//		std::cout << "Bad!";
	//	return;
	//}
	////polycr->OutputData();
	//if (ProcessMediumAngles(polycr))
	//{
	//	if (GlobalIntersectionCheck())
	//		std::cout << "Bad!";
	//	return;
	//}
	//polycr->OutputData();
	//ProcessLargeAngles(polycr);
	ProcessAngles(polycr);
	polycr->OutputData();
	if (GlobalIntersectionCheck())
		std::cout << "Bad!";
	SmoothMesh(10);
}

const bool Crystallite3::GlobalIntersectionCheck()
{
	for (auto &edge : innerEdges)
	{
		if (*edge &&
			EdgeGlobalIntersectionCheck(edge))
			return true;
	}

	return false;
}

void Crystallite3::SmoothMesh(int iterationsNum)
{
	for (int i = 0; i < iterationsNum; i++)
		for (auto &vert : innerVerts)
		{
			if (!*vert)
				continue;

			Vector3 shift;
			int delta_shifts_num = 0;
			for (auto &edge : innerEdges)
			{
				if (!*edge)
					continue;

				if (vert == (*edge)->vertexes[0])
				{
					shift += **(*edge)->vertexes[1] - **vert;
				}
				else if (vert == (*edge)->vertexes[1])
				{
					shift += **(*edge)->vertexes[0] - **vert;
				}
				delta_shifts_num++;
			}
			shift /= delta_shifts_num;
			(*vert)->SetPosition((*vert)->GetPosition() + shift);
		}
}

void Crystallite3::SmoothNotFinisedMesh(int iterationsNum)
{
	for (int i = 0; i < iterationsNum; i++)
		for (auto &vert : innerVerts)
		{
			if (!*vert)
				continue;

			SmoothAroundFrontVertex(vert);
		}
}

void Crystallite3::SmoothFront(int iterationsNum)
{
	for (int i = 0; i < iterationsNum; i++)
		for (auto &vert : innerVerts)
		{
			if (!*vert)
				continue;

			SmoothAroundFrontVertex(vert);
		}
}

void Crystallite3::SmoothAroundFrontVertex(unique_ptr<Vertex3>* frontVert)
{
	if ((*frontVert)->belongsToShellFacet ||
		(*frontVert)->belongsToShellEdge ||
		(*frontVert)->belongsToShellVertex)
	{
		return;
	}

	Vector3 shift;
	int delta_shifts_num = 0;
	for (auto &edge : innerEdges)
	{
		if (!*edge)
			continue;

		Vector3 d_shift;
		if (frontVert == (*edge)->vertexes[0])
		{
			d_shift = **(*edge)->vertexes[1] - **frontVert;
		}
		else if (frontVert == (*edge)->vertexes[1])
		{
			d_shift = **(*edge)->vertexes[0] - **frontVert;
		}
		shift += d_shift * (d_shift.Magnitude() - _preferredLength);
		delta_shifts_num++;
	}
	shift /= delta_shifts_num;
	(*frontVert)->SetPosition((*frontVert)->GetPosition() + shift);
}

void Crystallite3::AnalyzeMeshQuality(double &out_minQuality, double &out_avQuality)
{
	size_t simps_num = 0ull;
	double av_q = 0.0;
	double min_q = 1.0;
	for (auto &simp : innerSimps)
	{
		if (!*simp)
			continue;

		double q = (*simp)->Quality();
		av_q += q;
		simps_num++;
		if (q < min_q)
			min_q = q;
	}
	av_q /= simps_num;

	out_minQuality = min_q;
	out_avQuality = av_q;
}

Crystallite3::Crystallite3() {}

Crystallite3::~Crystallite3()
{
	for (auto &simp : innerSimps)
	{
		if (*simp)
			delete simp->release();
		delete simp;
	}
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