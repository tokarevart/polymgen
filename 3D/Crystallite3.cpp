#include "Crystallite3.h"
#include <algorithm>

#define DET(a, b, c, d) \
		(a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))


Crystallite3::Crystallite3() {}

Crystallite3::~Crystallite3()
{
	//for (auto &simp : innerSimps)
	//{
	//	delete simp->release();
	//	delete simp;
	//}
	for (auto &facet : innerFacets)
	{
		delete facet->release();
		delete facet;
	}
	for (auto &edge : innerEdges)
	{
		delete edge->release();
		delete edge;
	}
	for (auto &vert : innerVerts)
	{
		delete vert->release();
		delete vert;
	}
}