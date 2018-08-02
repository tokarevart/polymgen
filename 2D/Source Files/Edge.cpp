#include "Edge.h"

bool Edge::IsContaining(Node2 &node)
{
	if (nodes[0] == &node ||
		nodes[1] == &node)
	{
		return true;
	}

	return false;
}

void Edge::DestroyIfNoLinks()
{
	if (inclInSimplexes.empty() && !isShellPart)
	{
		delete this;
	}
}

Edge::Edge()
{
}

Edge::Edge(Node2 &node0, Node2 &node1)
{
	nodes[0] = &node0;
	nodes[1] = &node1;
}

Edge::~Edge()
{
	for (auto &nodei : nodes)
	{
		if (nodei)
		{
			nodei->inclInEdges.remove(this);
			for (auto nodej : nodes)
			{
				if (nodej != nodei)
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