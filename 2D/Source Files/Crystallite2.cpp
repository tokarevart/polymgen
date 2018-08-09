#include "Crystallite2.h"


const bool Crystallite2::IsContaining(const Simplex2& simp) const
{
	return find(simplexes.begin(), simplexes.end(), &simp) != simplexes.end();
}

const bool Crystallite2::IsContaining(const Edge2& edge) const
{
	for (auto simp : simplexes)
	{
		if (simp->IsContaining(edge))
		{
			return true;
		}
	}

	return false;
}

const bool Crystallite2::IsContaining(const Node2& node) const
{
	for (auto simp : simplexes)
	{
		if (simp->IsContaining(node))
		{
			return true;
		}
	}

	return false;
}

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
			delete *edge;
			delete edge;
			edge = nullptr;
		}
	}
	if (!simplexes.empty())
	{
		for (auto &simp : simplexes)
		{
			if (simp)
			{
				delete simp;
			}
		}
	}
}