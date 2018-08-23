#include "ShellEdge2.h"

const double MINUS_PI_DIV_2 = -3.141592653589793 * 0.5;

//void ShellEdge2::SetNormal()
//{
//	(*_normal)[0] = (*nodes[1]->get())[0] - (*nodes[0]->get())[0];
//	(*_normal)[1] = (*nodes[1]->get())[1] - (*nodes[0]->get())[1];
//
//	_normal->Rotate(MINUS_PI_DIV_2, Radian).Normalize();
//}

const double ShellEdge2::Magnitude() const
{
	return sqrt(SqrMagnitude());
}

const double ShellEdge2::SqrMagnitude() const
{
	Vector2 buf = **nodes[1] - **nodes[0];
	return Vector2::DotProduct(buf, buf);
}

void ShellEdge2::AttachNodes(const vector<unique_ptr<Node2>*>& free_nodes)
{
	Vector2 edge_direction = **nodes[1] - **nodes[0];

	unique_ptr<Node2>* nodes_attached_nodes[2];
	nodes_attached_nodes[0] = (*nodes[0])->FindAttachedNode(free_nodes);
	nodes_attached_nodes[1] = (*nodes[1])->FindAttachedNode(free_nodes);
	double max_cos = 0.0;
	unique_ptr<Node2>* node_with_max_cos = nullptr;
	for (auto &neighbor : (*nodes_attached_nodes[0])->neighbors)
	{
		if (neighbor == nodes_attached_nodes[1])
		{
			return;
		}
		double cosin = Vector2::Cos((*neighbor)->GetPosition(), (*nodes_attached_nodes[0])->GetPosition());
		if (cosin > max_cos)
		{
			max_cos = cosin;
			node_with_max_cos = neighbor;
		}
	}
	attachedNodes.push_back(node_with_max_cos);
	(*node_with_max_cos)->belongsToShellEdge = GetPtrToUniquePtr();
}

void ShellEdge2::ChangeAttachedNode(size_t index)
{
}

void ShellEdge2::SetAttachedNodesStartVectorsToEdge()
{
	size_t nodes_num = attachedNodes.size();
	attachedNodesStartVectorsToEdge.reserve(nodes_num);

	for (size_t i = 0; i < nodes_num; i++)
	{
		attachedNodesStartVectorsToEdge[i] =
			Vector2::Project(
				(*attachedNodes[i])->GetPosition(),
				(*nodes[0])->GetPosition(),
				(*nodes[1])->GetPosition())
			- (*attachedNodes[i])->GetPosition();
	}
}

void ShellEdge2::SetAttachedNodesDistanceFromStartPositionToEdge(const double& units, const double& outOf)
{
	double k = units / outOf;
	for (size_t i = 0, nodes_num = attachedNodes.size(); i < nodes_num; i++)
	{
		Vector2 cur_vec_to_edge = Vector2::Project(
		                               (*attachedNodes[i])->GetPosition(), 
		                               (*nodes[0])->GetPosition(), 
		                               (*nodes[1])->GetPosition());

		Vector2 delta_position = attachedNodesStartVectorsToEdge[i] * (k - 1.0) + cur_vec_to_edge;
		(*attachedNodes[i])->SetPosition((*attachedNodes[i])->GetPosition() + delta_position);
	}
}

const bool ShellEdge2::IsContaining(const ShellNode2& node) const
{
	if (nodes[0]->get() == &node ||
		nodes[1]->get() == &node)
	{
		return true;
	}

	return false;
}

ShellEdge2::ShellEdge2() : unique_ptr_helper<ShellEdge2>(this)
{
	nodes[0] = nullptr;
	nodes[1] = nullptr;
	//_normal.reset(new Vector2());
}

ShellEdge2::ShellEdge2(ShellNode2& node0, ShellNode2& node1) : unique_ptr_helper<ShellEdge2>(this)
{
	//_normal.reset(new Vector2());

	nodes[0] = node0.GetPtrToUniquePtr();
	nodes[1] = node1.GetPtrToUniquePtr();

	(*nodes[0])->inclInEdges.push_back(GetPtrToUniquePtr());
	(*nodes[1])->inclInEdges.push_back(GetPtrToUniquePtr());

	//SetNormal();
}

ShellEdge2::~ShellEdge2()
{
	for (auto &node : nodes)
	{
		if (*node)
		{
			if (std::find((*node)->inclInEdges.begin(), (*node)->inclInEdges.end(), GetPtrToUniquePtr()) != (*node)->inclInEdges.end())
			{
				(*node)->inclInEdges.remove(GetPtrToUniquePtr());
			}
			(*node)->DestroyIfNoLinks();
		}
	}
}