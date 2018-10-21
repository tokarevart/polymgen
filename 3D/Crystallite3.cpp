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

#define PI 3.141592653589793
#define PI_DIV_2 1.5707963267948966
#define DEG_80_IN_RADIANS 1.3962634015954636
#define DEG_90_IN_RADIANS PI_DIV_2
#define DEG_100_IN_RADIANS 1.7453292519943295
#define DEG_160_IN_RADIANS 2.7925268031909273
#define DEG_180_IN_RADIANS PI

template <class T>
void Crystallite3::ErasePtrsToNullptr(vector<unique_ptr<T>*> &vec)
{
	size_t real_objs_num = std::count_if(vec.begin(), vec.end(), [](unique_ptr<T>*& ptr) { return *ptr ? true : false; });
	vector<unique_ptr<T>*> buf_vec(real_objs_num);

	size_t firts_thr_nums = real_objs_num / 2ull;
	size_t second_thr_nums = real_objs_num - firts_thr_nums;
	//#pragma omp parallel num_threads(2)
	//{
	//#pragma omp single
	//{
	size_t index1 = 0ull;
	for (size_t i = 0ull; index1 < firts_thr_nums; i++)
		if (*vec[i])
		{
			buf_vec[index1] = vec[i];
			index1++;
		}
	//}
	//#pragma omp single
	//{
	size_t index2 = 0ull;
	for (size_t i = vec.size() - 1ull; index2 < second_thr_nums; i--)
		if (*vec[i])
		{
			buf_vec[real_objs_num - 1ull - index2] = vec[i];
			index2++;
		}
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

unique_ptr<FrontFacet3>* Crystallite3::FindFrontFacet(unique_ptr<Facet3>* &facet)
{
	for (auto &f_facet : frontFacets)
	{
		if (*f_facet &&
			(*f_facet)->facet == facet->get())
			return f_facet;
	}

	return nullptr;
}

unique_ptr<FrontEdge3>* Crystallite3::FindFrontEdge(const unique_ptr<Vertex3>* &v0, const unique_ptr<Vertex3>* &v1)
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

unique_ptr<FrontEdge3>* Crystallite3::FindFrontEdge(unique_ptr<Edge3>* &edge)
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
				+ (*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition())
			- v;

		for (auto &f_facet : frontFacets)
		{
			if (*f_facet && 
				Vector3::RayIntersectTriangle(
					v,
					dir,
					(*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition(),
					(*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition(),
					(*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition()))
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

const bool Crystallite3::EdgeIntersectionCheck(const unique_ptr<FrontEdge3>* frontEdge) // Needs improvement.
{
	unique_ptr<Vertex3>* opp_verts[2];
	(*frontEdge)->FindOppositeVertexes(
		frontFacets,
		frontEdges,
		opp_verts[0],
		opp_verts[1]);
	//double opp_poses_magn = ((*opp_verts[1])->GetPosition() - (*opp_verts[0])->GetPosition()).Magnitude();
	Vector3 opp_verts_poses[2];
	//opp_verts_poses[0] = (*opp_verts[0])->GetPosition() + Vector3(1.0, 1.4, 1.7) * 1e-4 * opp_poses_magn;
	//opp_verts_poses[1] = (*opp_verts[1])->GetPosition() - Vector3(1.4, 1.7, 1.0) * 1e-4 * opp_poses_magn;
	opp_verts_poses[0] = (*opp_verts[0])->GetPosition();
	opp_verts_poses[1] = (*opp_verts[1])->GetPosition();
	Vector3 opp_vert0_to_1 = opp_verts_poses[1] - opp_verts_poses[0];
	Vector3 delta_vec_0th_to_1st = 1e-3 * opp_vert0_to_1;

	if (LineSegmentFrontIntersectionCheck(
		opp_verts_poses[0] + delta_vec_0th_to_1st,
		opp_verts_poses[1] - delta_vec_0th_to_1st))
		return true;

	// Check distance between segments or point and segment.

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
					opp_verts_poses[1]) < 1e-5)
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
					opp_verts_poses[1]) < 1e-5)
				return true;
		}
		else
		{
			if (Vector3::SegmentsDistance(
					opp_verts_poses[0],
					opp_verts_poses[1],
					(*(*f_edge)->edge->vertexes[0])->GetPosition(),
					(*(*f_edge)->edge->vertexes[1])->GetPosition()) < 5e-3)
				return true;
		}
	}

	//for (auto &f_facet : frontFacets)
	//{
	//	if (!*f_facet ||
	//		!XOR((*f_facet)->facet->IsContaining(**opp_verts[0]), (*f_facet)->facet->IsContaining(**opp_verts[1])))
	//		continue;
	//
	//	Vector3 third = (*(*f_facet)->facet->FindVertexNotIncludedInEdge(**(*f_facet)->facet->edges[0]))->GetPosition();
	//	Vector3 normal = Vector3::CrossProduct(
	//		(*(*(*f_facet)->facet->edges[0])->vertexes[0])->GetPosition() - third,
	//		(*(*(*f_facet)->facet->edges[0])->vertexes[1])->GetPosition() - third);
	//
	//	if (abs(Vector3::DotProduct(normal, opp_vert0_to_1)) < 1e-8)
	//		return true;
	//}

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
	(*frontEdge)->FindFrontFacetsAround(
		frontFacets,
		f_facets_around[0], f_facets_around[1]);

	unique_ptr<Vertex3>* opp_verts[2];
	opp_verts[0] = (*f_facets_around[0])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);
	opp_verts[1] = (*f_facets_around[1])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);

	Vector3 plane0[2];
	plane0[0] = **opp_verts[0] - **(*frontEdge)->edge->vertexes[0];
	plane0[1] = **opp_verts[1] - **(*frontEdge)->edge->vertexes[0];
	Vector3 normal0 = Vector3::CrossProduct(plane0[0], plane0[1]);
	Vector3 plane1[2];
	plane1[0] = **opp_verts[0] - **(*frontEdge)->edge->vertexes[1];
	plane1[1] = **opp_verts[1] - **(*frontEdge)->edge->vertexes[1];
	Vector3 normal1 = Vector3::CrossProduct(plane1[0], plane1[1]);

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
			Vector3 f_normal = Vector3::CrossProduct(f_plane[0], f_plane[1]);

			if (abs(abs(Vector3::DotProduct(f_normal, normal0)) - sqrt(f_normal.SqrMagnitude() * normal0.SqrMagnitude())) < 1e-6 ||
				abs(abs(Vector3::DotProduct(f_normal, normal1)) - sqrt(f_normal.SqrMagnitude() * normal1.SqrMagnitude())) < 1e-6)
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

const bool Crystallite3::FrontContainsOfOnly1Simplex3OrEmpty()
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

void Crystallite3::CreateSimplex3AroundEdge(unique_ptr<FrontEdge3>* frontEdge)
{
	unique_ptr<FrontEdge3>* opp_f_edge = (*frontEdge)->FindOppositeFrontEdge(frontFacets, frontEdges);
	if (opp_f_edge)
	{
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

		unique_ptr<FrontEdge3>* new_f_facet_f_edges[3];
		unique_ptr<Edge3>* new_facet_edges[3];
		for (int i = 0; i < 3; i++)
		{
			new_facet_edges[i] = (*main_f_facets[i])->facet->FindEdgeNotContainingVertex(**main_vert);
			new_f_facet_f_edges[i] = FindFrontEdge(new_facet_edges[i]);
		}
		
		unique_ptr<FrontEdge3>* f_edge_ = nullptr;
		for (auto &f_facet : main_f_facets)
			for (auto &edge : (*f_facet)->facet->edges)
				if ((*edge)->IsContaining(**main_vert) &&
					(f_edge_ = FindFrontEdge(edge)))
					delete f_edge_->release();

		for (auto &f_facet : main_f_facets)
			delete f_facet->release();

		AddFrontFacet3(
			new_facet_edges[0],
			new_facet_edges[1],
			new_facet_edges[2]);

		// Maybe add new simplex3...
	}
	else
	{
		unique_ptr<FrontFacet3>* f_facets_around[2];
		(*frontEdge)->FindFrontFacetsAround(
			frontFacets,
			f_facets_around[0], f_facets_around[1]);

		unique_ptr<Vertex3>* opp_verts[2];
		opp_verts[0] = (*f_facets_around[0])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);
		opp_verts[1] = (*f_facets_around[1])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge);

		opp_f_edge = (new FrontEdge3(new Edge3(**opp_verts[0], **opp_verts[1])))->GetPtrToUniquePtr();
		AddFrontEdge3(opp_f_edge);

		unique_ptr<Edge3>* new_facet_edges[3];
		new_facet_edges[2] = (*opp_f_edge)->edge->GetPtrToUniquePtr();

		for (auto &edge : (*f_facets_around[0])->facet->edges)
			if ((*edge)->IsContaining(**(*frontEdge)->edge->vertexes[0]) &&
				((*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[0]) ||
				 (*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[1])))
			{
				new_facet_edges[0] = edge;
				break;
			}
		for (auto &edge : (*f_facets_around[1])->facet->edges)
			if ((*edge)->IsContaining(**(*frontEdge)->edge->vertexes[0]) &&
				((*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[0]) ||
				 (*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[1])))
			{
				new_facet_edges[1] = edge;
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
				break;
			}
		for (auto &edge : (*f_facets_around[1])->facet->edges)
			if ((*edge)->IsContaining(**(*frontEdge)->edge->vertexes[1]) &&
				((*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[0]) ||
				(*edge)->IsContaining(**(*opp_f_edge)->edge->vertexes[1])))
			{
				new_facet_edges[1] = edge;
				break;
			}
		AddFrontFacet3(
			new_facet_edges[0], 
			new_facet_edges[1], 
			new_facet_edges[2]);

		// Maybe add new simplex3...

		delete f_facets_around[0]->release();
		delete f_facets_around[1]->release();
		delete frontEdge->release();
	}
}

Vector3 Crystallite3::NewVertexPositionForSimplexesCreation(unique_ptr<FrontEdge3>* frontEdge, const bool calcVertDirByAngle)
{
	Vector3 f_edge_verts_poses[2];
	f_edge_verts_poses[0] = (*(*frontEdge)->edge->vertexes[0])->GetPosition();
	f_edge_verts_poses[1] = (*(*frontEdge)->edge->vertexes[1])->GetPosition();

	Vector3 orig = 0.5 * (f_edge_verts_poses[0] + f_edge_verts_poses[1]);

	unique_ptr<FrontFacet3>* f_facets[2];
	(*frontEdge)->FindFrontFacetsAround(frontFacets, f_facets[0], f_facets[1]);

	Vector3 f_facets_far_verts_poses[2];
	f_facets_far_verts_poses[0] = (*(*f_facets[0])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge))->GetPosition();
	f_facets_far_verts_poses[1] = (*(*f_facets[1])->facet->FindVertexNotIncludedInEdge(*(*frontEdge)->edge))->GetPosition();

	Vector3 facets_inner_vecs[2];
	facets_inner_vecs[0] = f_facets_far_verts_poses[0] - orig;
	facets_inner_vecs[1] = f_facets_far_verts_poses[1] - orig;

	Vector3 dir;
	if (calcVertDirByAngle &&
		(*frontEdge)->Angle(frontFacets) > PI)
		dir = -(facets_inner_vecs[0] + facets_inner_vecs[1]).Normalize();
	else
		dir = (facets_inner_vecs[0] + facets_inner_vecs[1]).Normalize();

	dir *= _preferredLength * 0.866; // sind(60) == 0.8660...
	Vector3 orig_plus_ex_dir = orig + dir * 2.0;
	while (
		LineSegmentFrontIntersectionCheck(f_edge_verts_poses[0] + 1e-3 * (orig_plus_ex_dir - f_edge_verts_poses[0]), orig_plus_ex_dir) ||
		LineSegmentFrontIntersectionCheck(f_edge_verts_poses[1] + 1e-3 * (orig_plus_ex_dir - f_edge_verts_poses[1]), orig_plus_ex_dir) ||
		LineSegmentFrontIntersectionCheck(f_facets_far_verts_poses[0] + 1e-3 * (orig_plus_ex_dir - f_facets_far_verts_poses[0]), orig_plus_ex_dir) ||
		LineSegmentFrontIntersectionCheck(f_facets_far_verts_poses[1] + 1e-3 * (orig_plus_ex_dir - f_facets_far_verts_poses[1]), orig_plus_ex_dir) ||
		!VertexInsideFrontCheck(orig_plus_ex_dir) || // Make facets with vertexes intersection check.
		//!VertexInsideFrontCheck(orig + dir) ||
		/*LineSegmentGlobalIntersectionCheck(f_edge_verts_poses[0] + 1e-3 * (orig + dir - f_edge_verts_poses[0]), orig + dir) ||
		LineSegmentGlobalIntersectionCheck(f_edge_verts_poses[1] + 1e-3 * (orig + dir - f_edge_verts_poses[1]), orig + dir) ||
		LineSegmentGlobalIntersectionCheck(f_facets_far_verts_poses[0] + 1e-3 * (orig + dir - f_facets_far_verts_poses[0]), orig + dir) ||
		LineSegmentFrontIntersectionCheck(f_facets_far_verts_poses[1] + 1e-3 * (orig + dir - f_facets_far_verts_poses[1]), orig + dir)*/)
	{
		dir *= 0.8;
		orig_plus_ex_dir = orig + dir * 2.0;
	}

	return orig + dir;
}

void Crystallite3::Create2Simplexes3AroundEdge(unique_ptr<FrontEdge3>* frontEdge, Vector3 vertPos)
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

	delete f_facets_around[0]->release();
	delete f_facets_around[1]->release();
	delete frontEdge->release();
}

template<class T>
size_t PtrsToNullptrNumber(vector<unique_ptr<T>*> &vec)
{
	return std::count_if(vec.begin(), vec.end(), [](unique_ptr<T>* &ptr) { return !*ptr ? true : false; });
}

const bool Crystallite3::ProcessSmallAngles()
{
	if (FrontContainsOfOnly1Simplex3OrEmpty())
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
		if (*frontEdges[i])
		{
			if ((*frontEdges[i])->Angle(frontFacets) < DEG_100_IN_RADIANS &&
				!FrontSplitCheck(frontEdges[i]) &&
				!EdgeIntersectionCheck(frontEdges[i]) &&
				!ParallelFacetsCheck(frontEdges[i]))
			{
				//std::cout << "\n\nAngle before creation: " << (*frontEdges[i])->Angle(frontFacets) << '\n';
				//std::cout << "Will be created:\n"
				//	<< "  " << Vector3((*(*frontEdges[i])->edge->vertexes[0])->GetPosition())[0] << '\n'
				//	<< "  " << Vector3((*(*frontEdges[i])->edge->vertexes[0])->GetPosition())[1] << '\n'
				//	<< "  " << Vector3((*(*frontEdges[i])->edge->vertexes[0])->GetPosition())[2] << '\n'
				//	<< "    " << Vector3((*(*frontEdges[i])->edge->vertexes[1])->GetPosition())[0] << '\n'
				//	<< "    " << Vector3((*(*frontEdges[i])->edge->vertexes[1])->GetPosition())[1] << '\n'
				//	<< "    " << Vector3((*(*frontEdges[i])->edge->vertexes[1])->GetPosition())[2] << '\n';
				CreateSimplex3AroundEdge(frontEdges[i]);
				//polycr->OutputData();
				i = 0ull; //return;
				
				//int copys = 0;
				//for (auto &edgei : innerEdges)
				//	for (auto &edgej : innerEdges)
				//	{
				//		if (edgei != edgej &&
				//			((*edgei)->vertexes[0] == (*edgej)->vertexes[0] &&
				//			(*edgei)->vertexes[1] == (*edgej)->vertexes[1]) ||
				//			((*edgei)->vertexes[0] == (*edgej)->vertexes[1] &&
				//			(*edgei)->vertexes[1] == (*edgej)->vertexes[0]))
				//		{
				//			copys++;
				//			//std::cout << "\nEdge copy:\n"
				//			//	<< "  " << Vector3((*(*edgei)->vertexes[0])->GetPosition())[0] << '\n'
				//			//	<< "  " << Vector3((*(*edgei)->vertexes[0])->GetPosition())[1] << '\n'
				//			//	<< "  " << Vector3((*(*edgei)->vertexes[0])->GetPosition())[2] << '\n'
				//			//	<< "    " << Vector3((*(*edgei)->vertexes[1])->GetPosition())[0] << '\n'
				//			//	<< "    " << Vector3((*(*edgei)->vertexes[1])->GetPosition())[1] << '\n'
				//			//	<< "    " << Vector3((*(*edgei)->vertexes[1])->GetPosition())[2] << '\n'
				//			//	<< " " << "Angle: " << (*frontEdges[i])->Angle(frontFacets) << '\n';
				//		}
				//	}
				//if (copys > 1)
				//	std::cout << "Copys: " << copys;

				if (FrontContainsOfOnly1Simplex3OrEmpty())
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
				//return;

				if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 3ull)
					ErasePtrsToNullptr(frontEdges);
				if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 3ull)
					ErasePtrsToNullptr(frontFacets);
			}
			else i++;
		}
		else i++;
	}

	return false;
}

void Crystallite3::ProcessMediumAngles(Polycrystalline3* polycr)
{
	if (FrontContainsOfOnly1Simplex3OrEmpty())
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
		if (!*frontEdges[i] ||
			FrontSplitCheck(frontEdges[i]) ||
			EdgeIntersectionCheck(frontEdges[i]) ||
			ParallelFacetsCheck(frontEdges[i]) ||
			(*frontEdges[i])->Angle(frontFacets) >= DEG_180_IN_RADIANS)
		{
			i++;
			continue;
		}

		Create2Simplexes3AroundEdge(frontEdges[i], NewVertexPositionForSimplexesCreation(frontEdges[i]));
		polycr->OutputData();
		if (ProcessSmallAngles())
			return;
		polycr->OutputData();
		i = 0ull;

		if (FrontContainsOfOnly1Simplex3OrEmpty())
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

		if (PtrsToNullptrNumber(frontEdges) > frontEdges.size() / 3ull)
			ErasePtrsToNullptr(frontEdges);
		if (PtrsToNullptrNumber(frontFacets) > frontFacets.size() / 3ull)
			ErasePtrsToNullptr(frontFacets);
	}

	polycr->OutputData();
	if (!FrontContainsOfOnly1Simplex3OrEmpty())
		throw std::logic_error("Crystallite volume wasn't exhausted due to unexpected logical error.");
}

void Crystallite3::TriangulateVolume(const double preferredLength, Polycrystalline3* polycr)
{
	_preferredLength = preferredLength;
	polycr->OutputData();
	if (ProcessSmallAngles())
		return;
	polycr->OutputData();
	ProcessMediumAngles(polycr);
	polycr->OutputData();
	if (GlobalIntersectionCheck())
		std::cout << "Bad!";
	//ProcessLargeAngles();
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