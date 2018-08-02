#include "Crystallite2.h"


bool Crystallite2::IsContaining(Simplex2 &simp)
{
	return find(simplexes.begin(), simplexes.end(), &simp) != simplexes.end();
}

bool Crystallite2::IsContaining(Edge &edge)
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

bool Crystallite2::IsContaining(Node2 &node)
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

bool Crystallite2::Contains(Simplex2 &simp)
{
	// Do something

	return false;
}

bool Crystallite2::Contains(Edge &edge)
{
	// Do something

	return false;
}

bool Crystallite2::Contains(Node2 &node)
{
	// Do something

	return false;
}

Crystallite2::Crystallite2()
{
}

Crystallite2::~Crystallite2()
{
	for (auto &edge : shellEdges)
	{
		if (edge)
		{
			delete edge;
		}
	}
	for (auto &simp : simplexes)
	{
		if (simp)
		{
			delete simp;
		}
	}
}