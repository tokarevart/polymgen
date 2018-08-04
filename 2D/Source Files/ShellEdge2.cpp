#include "ShellEdge2.h"

const double MINUS_PI_DIV_2 = -3.141592653589793 * 0.5;

void ShellEdge2::SetNormal()
{
	(*_normal)[0] = (*nodes[1])[0] - (*nodes[0])[0];
	(*_normal)[1] = (*nodes[1])[1] - (*nodes[0])[1];

	_normal->Rotate(MINUS_PI_DIV_2, Radian).Normalize();
}

const bool ShellEdge2::IsContaining(const ShellNode2& node)
{
	if (nodes[0] == &node ||
		nodes[1] == &node)
	{
		return true;
	}

	return false;
}

const bool ShellEdge2::IsAttached(const Node2& node)
{
	return find(attachedNodes.begin(), attachedNodes.end(), &node) != attachedNodes.end();
}

ShellEdge2::ShellEdge2()
{
}

ShellEdge2::ShellEdge2(ShellNode2& node0, ShellNode2& node1)
{
	nodes[0] = &node0;
	nodes[1] = &node1;

	nodes[0]->inclInEdges.push_back(this);
	nodes[1]->inclInEdges.push_back(this);

	SetNormal();
}

ShellEdge2::~ShellEdge2()
{
	for (auto &node : nodes)
	{
		if (node)
		{
			node->inclInEdges.remove(this);
		}
	}
	for (auto attedNode : attachedNodes)
	{
		if (attedNode)
		{
			attedNode->belongsToShellEdge = nullptr;
		}
	}
	delete _normal;
}