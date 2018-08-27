#include "Crystallite2.h"
#include <algorithm>

#define DET(a, b, c, d) \
		(a * d - b * c)

#define EPS 1e-10
#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))


const bool Crystallite2::Contains(const Node2& node) const
{
	int hits_num[3] = { 0, 0, 0 };
	Vector2 intersect_point;
	for (auto &s_edge : shellEdges)
	{
		Vector2 ray0(1.0, 0.2);
		intersect_point =
			Vector2::LinesIntersection(
				node.GetPosition(),
				node.GetPosition() + ray0,
				(*(*s_edge)->nodes[0])->GetPosition(),
				(*(*s_edge)->nodes[1])->GetPosition());
		if (INSIDE_RECTANGLE(
			(*(*s_edge)->nodes[0])->GetPosition(),
			(*(*s_edge)->nodes[1])->GetPosition(),
			intersect_point) &&
			Vector2::DotProduct(intersect_point - node.GetPosition(), ray0) > 0.0)
		{
			hits_num[0]++;
		}


		Vector2 ray1(0.2, -1.0);
		intersect_point =
			Vector2::LinesIntersection(
				node.GetPosition(),
				node.GetPosition() + ray1,
				(*(*s_edge)->nodes[0])->GetPosition(),
				(*(*s_edge)->nodes[1])->GetPosition());
		if (INSIDE_RECTANGLE(
			(*(*s_edge)->nodes[0])->GetPosition(),
			(*(*s_edge)->nodes[1])->GetPosition(),
			intersect_point) &&
			Vector2::DotProduct(intersect_point - node.GetPosition(), ray1) > 0.0)
		{
			hits_num[1]++;
		}


		Vector2 ray2(-1.0, 1.0);
		intersect_point =
			Vector2::LinesIntersection(
				node.GetPosition(),
				node.GetPosition() + ray2,
				(*(*s_edge)->nodes[0])->GetPosition(),
				(*(*s_edge)->nodes[1])->GetPosition());
		if (INSIDE_RECTANGLE(
			(*(*s_edge)->nodes[0])->GetPosition(),
			(*(*s_edge)->nodes[1])->GetPosition(),
			intersect_point) &&
			Vector2::DotProduct(intersect_point - node.GetPosition(), ray2) > 0.0)
		{
			hits_num[2]++;
		}
	}

	hits_num[0] %= 2;
	hits_num[1] %= 2;
	hits_num[2] %= 2;
	if (hits_num[0] == hits_num[1] &&
		hits_num[1] == hits_num[2] &&
		hits_num[2] == hits_num[0])
	{
		if (hits_num[0] == 1)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	else
	{
		if ((hits_num[0] == hits_num[1] && hits_num[1] == 1) ||
			(hits_num[0] == hits_num[2] && hits_num[2] == 1) ||
			(hits_num[1] == hits_num[2] && hits_num[2] == 1))
		{
			return true;
		}
		else
		{
			return false;
		}
	}
}

Crystallite2::Crystallite2() {}

Crystallite2::~Crystallite2()
{
	for (auto &edge : shellEdges)
	{
		if (*edge)
		{
			delete edge->release();
		}
	}
}