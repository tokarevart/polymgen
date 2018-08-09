#include "Edge2.h"

const double Edge2::Magnitude() const
{
	Vector2 buf = **nodes[1] - **nodes[0];
	return buf.Magnitude();
}

const double Edge2::SqrMagnitude() const
{
	Vector2 buf = **nodes[1] - **nodes[0];
	return buf.SqrMagnitude();
}

const bool Edge2::IsContaining(const Node2& node)
{
	if (nodes[0]->get() == &node ||
		nodes[1]->get() == &node)
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

Edge2::Edge2() : unique_ptr_helper(this) {}

Edge2::Edge2(Node2& node0, Node2& node1) : unique_ptr_helper(this)
{
	nodes[0] = node0.GetPtrToUniquePtr();
	nodes[1] = node1.GetPtrToUniquePtr();

	(*nodes[0])->inclInEdges.push_back(GetPtrToUniquePtr());
	(*nodes[0])->neighbors.push_back(nodes[1]);

	(*nodes[1])->inclInEdges.push_back(GetPtrToUniquePtr());
	(*nodes[1])->neighbors.push_back(nodes[0]);
}

Edge2::~Edge2()
{
	for (auto &nodei : nodes)
	{
		if (*nodei)
		{
			(*nodei)->inclInEdges.remove(GetPtrToUniquePtr());
			for (auto &nodej : nodes)
			{
				if (*nodej && (nodej != nodei))
				{
					(*nodei)->neighbors.remove(nodej);
				}
			}
			(*nodei)->DestroyIfNoLinks();
		}
	}

	for (auto &simp : inclInSimplexes)
	{
		if (*simp)
		{
			simp->release();
		}
	}
}