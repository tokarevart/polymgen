#include "Crystallite3.h"
#include <algorithm>

#define DET(a, b, c, d) \
		(a * d - b * c)

#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) < p && p < std::max(p0_coor, p1_coor))

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))


const bool Crystallite3::Contains(const Vertex3& vertex) const
{
}

Crystallite3::Crystallite3() {}

Crystallite3::~Crystallite3()
{
}