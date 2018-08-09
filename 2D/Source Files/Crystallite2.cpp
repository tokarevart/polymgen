#include "Crystallite2.h"


const bool Crystallite2::Contains(const Simplex2& simp) const
{
	// Do something

	return false;
}

const bool Crystallite2::Contains(const Edge2& edge) const
{
	// Do something

	return false;
}

const bool Crystallite2::Contains(const Node2& node) const
{
	// Do something

	return false;
}

Crystallite2::Crystallite2() {}

Crystallite2::~Crystallite2()
{
	for (auto &edge : shellEdges)
	{
		if (*edge)
		{
			edge->release();
			delete edge;
		}
	}
	for (auto &simp : simplexes)
	{
		if (*simp)
		{
			simp->release();
			delete simp;
		}
	}
}