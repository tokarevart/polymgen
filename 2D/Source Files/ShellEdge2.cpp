#include "ShellEdge2.h"

#define PI_DIV_2 1.5707963267948966
#define MINUS_PI_DIV_2 -1.5707963267948966


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
	unique_ptr<Node2>* current_node = nodes_attached_nodes[0];
	for (size_t i = 0, max = free_nodes.size(); i < max; i++)
	{
		for (auto &neighbor : (*current_node)->neighbors)
		{
			if (neighbor == nodes_attached_nodes[1])
			{
				return;
			}
			double cos_alpha = Vector2::Cos(**neighbor - **nodes_attached_nodes[0], edge_direction);
			double beta = acos(Vector2::Cos(**neighbor - **current_node, **nodes_attached_nodes[1] - **current_node));
			double gamma = acos(Vector2::Cos(**nodes_attached_nodes[1] - **current_node, edge_direction));
			if (beta + gamma < PI_DIV_2 &&
				cos_alpha > max_cos)
			{
				max_cos = cos_alpha;
				node_with_max_cos = neighbor;
			}
		}
		attachedNodes.push_back(node_with_max_cos);
		(*node_with_max_cos)->belongsToShellEdge = GetPtrToUniquePtr();
		current_node = node_with_max_cos;
		max_cos = 0.0;
		node_with_max_cos = nullptr;
	}

	throw std::exception("Something went wrong in ShellEdge2::AttachNodes");
}

void ShellEdge2::ChangeAttachedNode(const size_t& index, const vector<unique_ptr<Node2>*>& free_nodes)
{
	unique_ptr<Node2>* new_node = nullptr;
	unique_ptr<Node2>* shell_node = nullptr;
	unique_ptr<Node2>* around_node = nullptr;
	if (index < attachedNodes.size() / 2)
	{
		shell_node = (*nodes[0])->FindAttachedNode(free_nodes);
		around_node = attachedNodes[index + 1];
	}
	else 
	{
		shell_node = (*nodes[1])->FindAttachedNode(free_nodes);
		if (index > 0)
		{
			around_node = attachedNodes[index - 1];
		}
		else
		{
			around_node = (*nodes[0])->FindAttachedNode(free_nodes);
		}
	}
	
	for (auto &shell_node_neighbor : (*shell_node)->neighbors)
	{
		for (auto &around_node_neighbor : (*around_node)->neighbors)
		{
			if (shell_node_neighbor == around_node_neighbor &&
				shell_node_neighbor != attachedNodes[index])
			{
				attachedNodes[index] = shell_node_neighbor;
				(*attachedNodes[index])->belongsToShellEdge = GetPtrToUniquePtr();
				return;
			}
		}
	}
}

void ShellEdge2::SetAttachedNodesStartVectorsToEdge()
{
	size_t nodes_num = attachedNodes.size();
	if (nodes_num > 0)
	{
		attachedNodesStartVectorsToEdge.assign(nodes_num, Vector2());
	}

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
		Vector2 cur_vec_to_edge = 
			Vector2::Project(
				(*attachedNodes[i])->GetPosition(), 
				(*nodes[0])->GetPosition(), 
				(*nodes[1])->GetPosition())
			- (*attachedNodes[i])->GetPosition();

		Vector2 delta_position = attachedNodesStartVectorsToEdge[i] * (k - 1.0) + cur_vec_to_edge;
		(*attachedNodes[i])->SetPosition((*attachedNodes[i])->GetPosition() + delta_position);
	}
}

const bool ShellEdge2::ContainsAttachedNode(const unique_ptr<Node2>* const& node)
{
	for (auto &attached_node : attachedNodes)
	{
		if (attached_node == node)
		{
			return true;
		}
	}

	return false;
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
}

ShellEdge2::ShellEdge2(ShellNode2& node0, ShellNode2& node1) : unique_ptr_helper<ShellEdge2>(this)
{
	nodes[0] = node0.GetPtrToUniquePtr();
	nodes[1] = node1.GetPtrToUniquePtr();

	(*nodes[0])->inclInEdges.push_back(GetPtrToUniquePtr());
	(*nodes[1])->inclInEdges.push_back(GetPtrToUniquePtr());
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