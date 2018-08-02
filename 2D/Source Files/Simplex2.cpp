#include "Simplex2.h"


bool Simplex2::IsContaining(Edge &edge)
{
	for (int i = 0; i < 3; i++)
	{
		if (edges[i] == &edge)
		{
			return true;
		}
	}

	return false;
}

bool Simplex2::IsContaining(Node2 &node)
{
	for (auto &edge : edges)
	{
		if (edge->IsContaining(node))
		{
			return true;
		}
	}

	return false;
}

Simplex2::Simplex2()
{
}

Simplex2::Simplex2(Edge edge0, Edge edge1, Edge edge2)
{
	edges[0] = &edge0;
	edges[1] = &edge1;
	edges[2] = &edge2;
}

Simplex2::~Simplex2()
{
	for (auto &edge : edges)
	{
		if (edge)
		{
			edge->inclInSimplexes.remove(this);
			edge->DestroyIfNoLinks();
		}
	}
}