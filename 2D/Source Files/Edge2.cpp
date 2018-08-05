#include "Edge2.h"

const double Edge2::GetLength()
{
	Vector2 buf = (*nodes)[1] - (*nodes)[0];
	return buf.GetLength();
}

const double Edge2::GetSqrLength()
{
	Vector2 buf = (*nodes)[1] - (*nodes)[0];
	return buf.GetSqrLength();
}

const bool Edge2::IsContaining(const Node2 &node)
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

	nodes[0]->inclInEdges.push_back(this);
	nodes[0]->neighbors.push_back(nodes[1]);

	nodes[1]->inclInEdges.push_back(this);
	nodes[1]->neighbors.push_back(nodes[0]);
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