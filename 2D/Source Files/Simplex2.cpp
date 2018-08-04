#include "Simplex2.h"


bool Simplex2::IsContaining(Edge2 &edge)
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

Simplex2::Simplex2(Edge2 &edge0, Edge2 &edge1, Edge2 &edge2)
{
	edges[0] = &edge0;
	edges[1] = &edge1;
	edges[2] = &edge2;
}

Simplex2::Simplex2(Node2 &node0, Node2 &node1, Node2 &node2)
{
	edges[0] = new Edge2(node0, node1);
	edges[1] = new Edge2(node1, node2);
	edges[2] = new Edge2(node2, node0);
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