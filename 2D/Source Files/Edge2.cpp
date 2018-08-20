#include "Edge2.h"
#include <algorithm>


#define EPS 1e-10
#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))


const double Edge2::Magnitude() const
{
	return (**nodes[1] - **nodes[0]).Magnitude();
}

const double Edge2::SqrMagnitude() const
{
	return (**nodes[1] - **nodes[0]).SqrMagnitude();
}

unique_ptr<Node2>* FindSimplexNodeNotBelongsToEdge(Simplex2& simp, Edge2& edge)
{
	for (auto &simp_edge : simp.edges)
	{
		if (simp_edge->get() != &edge)
		{
			for (auto &node : (*simp_edge)->nodes)
			{
				if (!edge.IsContaining(**node))
				{
					return node;
				}
			}
		}
	}

	return nullptr;
}

void Edge2::MakeTwoInstead(list<unique_ptr<Simplex2>*>& freeSimplexes, vector<unique_ptr<Edge2>*>& freeEdges, vector<unique_ptr<Node2>*>& freeNodes)
{
	unique_ptr<Node2>* inner_node = (new Node2((*nodes[0])->GetPosition() + 0.5 * (**nodes[1] - **nodes[0])))->GetPtrToUniquePtr();
	if ((*nodes[0]))
	{
		list<Crystallite2*> inner_node_cryses;
		for (auto &crys : (*nodes[0])->belongsToCryses)
		{
			if (std::find(
					(*nodes[1])->belongsToCryses.begin(),
					(*nodes[1])->belongsToCryses.end(),
					crys)
				!= (*nodes[1])->belongsToCryses.end())
			{
				inner_node_cryses.push_back(crys);
			}
		}
		(*inner_node)->belongsToCryses = inner_node_cryses;
	}
	if (BelongsToShell())
	{
		for (auto &node : nodes)
		{
			if (*node && (*node)->belongsToShellEdge)
			{
				(*inner_node)->belongsToShellEdge = (*node)->belongsToShellEdge;
			}
		}
	}
	freeNodes.push_back(inner_node);

	vector<unique_ptr<Simplex2>*> oldSimps;
	for (auto &simp : inclInSimplexes)
	{
		if (*simp)
		{
			oldSimps.push_back(simp);
		}
	}

	unique_ptr<Edge2>* edge_halfs[2];
	edge_halfs[0] = (new Edge2(**nodes[0], **inner_node))->GetPtrToUniquePtr();
	edge_halfs[1] = (new Edge2(**inner_node, **nodes[1]))->GetPtrToUniquePtr();

	freeEdges.push_back(edge_halfs[0]);
	freeEdges.push_back(edge_halfs[1]);

	unique_ptr<Node2>* not_belongs_to_edge_node;
	unique_ptr<Edge2>* dividing_edge;
	for (auto &simp : oldSimps)
	{
		not_belongs_to_edge_node = FindSimplexNodeNotBelongsToEdge(**simp, *this);
		dividing_edge = (new Edge2(**inner_node, **not_belongs_to_edge_node))->GetPtrToUniquePtr();
		freeEdges.push_back(dividing_edge);

		freeSimplexes.push_back(
			(new Simplex2(
				**(*simp)->FindEdge(**nodes[0], **not_belongs_to_edge_node),
				**edge_halfs[0],
				**dividing_edge))
			->GetPtrToUniquePtr());

		freeSimplexes.push_back(
			(new Simplex2(
				**(*simp)->FindEdge(**nodes[1], **not_belongs_to_edge_node),
				**edge_halfs[1],
				**dividing_edge))
			->GetPtrToUniquePtr());

		delete simp->release();
	}
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

const bool Edge2::BelongsToShell()
{
	if (((*nodes[0])->belongsToShellNode && (*nodes[1])->belongsToShellNode) ||
		((*nodes[0])->belongsToShellEdge && (*nodes[1])->belongsToShellNode) ||
		((*nodes[0])->belongsToShellEdge == (*nodes[1])->belongsToShellEdge && (*nodes[1])->belongsToShellEdge))
	{
		return true;
	}

	return false;
}

const bool Edge2::IntersectsWith(const Edge2& edge)
{
	if (Vector2::SegmentsIntersection(
			(*nodes[0])->GetPosition(), 
			(*nodes[1])->GetPosition(), 
			(*edge.nodes[0])->GetPosition(), 
			(*edge.nodes[1])->GetPosition()))
	{
		return true;
	}

	return false;
}

const bool Edge2::IntersectsWith(const ShellEdge2& edge)
{
	if (Vector2::SegmentsIntersection(
			(*nodes[0])->GetPosition(),
			(*nodes[1])->GetPosition(),
			(*edge.nodes[0])->GetPosition(),
			(*edge.nodes[1])->GetPosition()))
	{
		return true;
	}

	return false;
}

const bool Edge2::IntersectsWith(const vector<unique_ptr<ShellEdge2>*>& shellEdges)
{
	for (auto &s_edge : shellEdges)
	{
		if (*s_edge && IntersectsWith(**s_edge))
		{
			return true;
		}
	}

	return false;
}

void Edge2::DestroyIfNoLinks()
{
	if (inclInSimplexes.empty())
	{
		delete _uniquePtr->release();
	}
}

Edge2::Edge2() : unique_ptr_helper<Edge2>(this) 
{
	for (auto &node : nodes)
	{
		node = nullptr;
	}
}

Edge2::Edge2(Node2& node0, Node2& node1) : unique_ptr_helper<Edge2>(this)
{
	nodes[0] = node0.GetPtrToUniquePtr();
	nodes[1] = node1.GetPtrToUniquePtr();

	(*nodes[0])->inclInEdges.push_back(GetPtrToUniquePtr());
	if (std::find(
			(*nodes[0])->neighbors.begin(), 
			(*nodes[0])->neighbors.end(), 
			nodes[1]) 
		== (*nodes[0])->neighbors.end())
	{
		(*nodes[0])->neighbors.push_back(nodes[1]);
	}

	(*nodes[1])->inclInEdges.push_back(GetPtrToUniquePtr());
	if (std::find(
			(*nodes[1])->neighbors.begin(), 
			(*nodes[1])->neighbors.end(), 
			nodes[0]) 
		== (*nodes[1])->neighbors.end())
	{
		(*nodes[1])->neighbors.push_back(nodes[0]);
	}
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
			delete simp->release();
		}
	}
}