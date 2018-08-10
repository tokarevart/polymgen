#include "Crystallite2.h"


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
			delete edge->release();
		}
	}

	for (auto &simp : simplexes)
	{
		if (*simp)
		{
			delete simp->release();
		}
	}
}