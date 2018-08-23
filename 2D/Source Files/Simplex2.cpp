#include "Simplex2.h"


//Simplex2& Simplex2::UpdateAverageEdgesLength()
//{
//	averageEdgesLength = ((*edges[0])->Magnitude() + (*edges[1])->Magnitude() + (*edges[2])->Magnitude()) * 0.3333333333333333;
//
//	return *this;
//}

unique_ptr<Edge2>* Simplex2::FindEdge(Node2& node0, Node2& node1)
{
	for (auto &edge : edges)
	{
		if (((*edge)->nodes[0]->get() == &node0 &&
			 (*edge)->nodes[1]->get() == &node1) ||
			((*edge)->nodes[0]->get() == &node1 &&
			 (*edge)->nodes[1]->get() == &node0))
		{
			return edge;
		}
	}
	
	return nullptr;
}

unique_ptr<Edge2>* Simplex2::MinEdge()
{
	if ((*edges[0])->SqrMagnitude() < (*edges[1])->SqrMagnitude())
	{
		if ((*edges[0])->SqrMagnitude() < (*edges[2])->SqrMagnitude())
		{
			return edges[0];
		}
		else
		{
			return edges[2];
		}
	}
	else
	{
		if ((*edges[1])->SqrMagnitude() < (*edges[2])->SqrMagnitude())
		{
			return edges[1];
		}
		else
		{
			return edges[2];
		}
	}

	return nullptr;
}

unique_ptr<Edge2>* Simplex2::MaxEdge()
{
	if ((*edges[0])->SqrMagnitude() > (*edges[1])->SqrMagnitude())
	{
		if ((*edges[0])->SqrMagnitude() > (*edges[2])->SqrMagnitude())
		{
			return edges[0];
		}
		else
		{
			return edges[2];
		}
	}
	else
	{
		if ((*edges[1])->SqrMagnitude() > (*edges[2])->SqrMagnitude())
		{
			return edges[1];
		}
		else
		{
			return edges[2];
		}
	}

	return nullptr;
}

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

Simplex2::Simplex2() : unique_ptr_helper<Simplex2>(this) 
{
	for (auto &edge : edges)
	{
		edge = nullptr;
	}
}

Simplex2::Simplex2(Edge2& edge0, Edge2& edge1, Edge2& edge2) : unique_ptr_helper<Simplex2>(this)
{
	edges[0] = edge0.GetPtrToUniquePtr();
	edges[1] = edge1.GetPtrToUniquePtr();
	edges[2] = edge2.GetPtrToUniquePtr();

	//averageEdgesLength = (edge0.Magnitude() + edge1.Magnitude() + edge2.Magnitude()) * 0.3333333333333333;

	(*edges[0])->inclInSimplexes.push_back(GetPtrToUniquePtr());
	(*edges[1])->inclInSimplexes.push_back(GetPtrToUniquePtr());
	(*edges[2])->inclInSimplexes.push_back(GetPtrToUniquePtr());

	for (auto &edge : edges)
	{
		for (auto &node : (*edge)->nodes)
		{
			if (std::find(
					(*node)->inclInSimplexes.begin(),
					(*node)->inclInSimplexes.end(), 
					GetPtrToUniquePtr())
				== (*node)->inclInSimplexes.end())
			{
				(*node)->inclInSimplexes.push_back(GetPtrToUniquePtr());
			}
		}
	}
}

Simplex2::Simplex2(Node2& node0, Node2& node1, Node2& node2) : unique_ptr_helper<Simplex2>(this)
{
	edges[0] = (new Edge2(node0, node1))->GetPtrToUniquePtr();
	edges[1] = (new Edge2(node1, node2))->GetPtrToUniquePtr();
	edges[2] = (new Edge2(node2, node0))->GetPtrToUniquePtr();

	//averageEdgesLength = ((*edges[0])->Magnitude() + (*edges[1])->Magnitude() + (*edges[2])->Magnitude()) * 0.3333333333333333;

	(*edges[0])->inclInSimplexes.push_back(GetPtrToUniquePtr());
	(*edges[1])->inclInSimplexes.push_back(GetPtrToUniquePtr());
	(*edges[2])->inclInSimplexes.push_back(GetPtrToUniquePtr());

	for (auto &edge : edges)
	{
		for (auto &node : (*edge)->nodes)
		{
			if (std::find(
					(*node)->inclInSimplexes.begin(),
					(*node)->inclInSimplexes.end(),
					GetPtrToUniquePtr())
				== (*node)->inclInSimplexes.end())
			{
				(*node)->inclInSimplexes.push_back(GetPtrToUniquePtr());
			}
		}
	}
}

Simplex2::~Simplex2()
{
	for (auto &edge : edges)
	{
		if (*edge)
		{
			for (auto &node : (*edge)->nodes)
			{
				if (*node)
				{
					(*node)->inclInSimplexes.remove(GetPtrToUniquePtr());
				}
			}

			(*edge)->inclInSimplexes.remove(GetPtrToUniquePtr());
			(*edge)->DestroyIfNoLinks();
		}
	}
}