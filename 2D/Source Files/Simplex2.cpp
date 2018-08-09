#include "Simplex2.h"


const bool Simplex2::IsContaining(const Edge2& edge) const
{
	for (int i = 0; i < 3; i++)
	{
		if (edges[i]->get() == &edge)
		{
			return true;
		}
	}

	return false;
}

const bool Simplex2::IsContaining(const Node2& node) const
{
	for (auto &edge : edges)
	{
		if ((*edge)->IsContaining(node))
		{
			return true;
		}
	}

	return false;
}

Simplex2::Simplex2() : unique_ptr_helper(this) {}

Simplex2::Simplex2(Edge2& edge0, Edge2& edge1, Edge2& edge2) : unique_ptr_helper(this)
{
	edges[0] = edge0.GetPtrToUniquePtr();
	edges[1] = edge1.GetPtrToUniquePtr();
	edges[2] = edge2.GetPtrToUniquePtr();

	(*edges[0])->inclInSimplexes.push_back(GetPtrToUniquePtr());
	(*edges[1])->inclInSimplexes.push_back(GetPtrToUniquePtr());
	(*edges[2])->inclInSimplexes.push_back(GetPtrToUniquePtr());
}

Simplex2::Simplex2(Node2& node0, Node2& node1, Node2& node2) : unique_ptr_helper(this)
{
	edges[0] = new Edge2(node0, node1)->GetPtrToUniquePtr();
	edges[1] = new Edge2(node1, node2)->GetPtrToUniquePtr();
	edges[2] = new Edge2(node2, node0)->GetPtrToUniquePtr();

	(*edges[0])->inclInSimplexes.push_back(GetPtrToUniquePtr());
	(*edges[1])->inclInSimplexes.push_back(GetPtrToUniquePtr());
	(*edges[2])->inclInSimplexes.push_back(GetPtrToUniquePtr());
}

Simplex2::~Simplex2()
{
	for (auto &edge : edges)
	{
		if (*edge)
		{
			(*edge)->inclInSimplexes.remove(GetPtrToUniquePtr());
			(*edge)->DestroyIfNoLinks();
		}
	}
}