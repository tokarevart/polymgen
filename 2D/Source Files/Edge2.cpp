#include "Edge2.h"

bool Edge2::IsContaining(Node2 &node)
{
	if (nodes[0] == &node ||
		nodes[1] == &node)
	{
		return true;
	}

	return false;
}

void Edge2::DestroyIfNoLinks()
{
	if (inclInSimplexes.empty())
	{
		delete this;
	}
}

Edge2::Edge2()
{
}

Edge2::Edge2(Node2 &node0, Node2 &node1)
{
	nodes[0] = &node0;
	nodes[1] = &node1;
}

Edge2::~Edge2()
{
	for (auto &nodei : nodes)
	{
		if (nodei)
		{
			nodei->inclInEdges.remove(this);
			for (auto nodej : nodes)
			{
				if (nodej && nodej != nodei)
				{
					nodei->neighbors.remove(nodej);
				}
			}
			nodei->DestroyIfNoLinks();
		}
	}

	for (auto &simp : inclInSimplexes)
	{
		if (simp)
		{
			delete simp;
		}
	}
}