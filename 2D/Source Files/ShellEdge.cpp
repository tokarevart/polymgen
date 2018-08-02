#include "ShellEdge.h"


bool ShellEdge::Contains(ShellNode2 &node)
{
	if (nodes[0] == &node ||
		nodes[1] == &node)
	{
		return true;
	}

	return false;
}

bool ShellEdge::IsAttached(Node2 &node)
{
	return find(attachedNodes.begin(), attachedNodes.end(), &node) != attachedNodes.end();
}

ShellEdge::ShellEdge()
{
}

ShellEdge::ShellEdge(ShellNode2 &node0, ShellNode2 &node1)
{
	nodes[0] = &node0;
	nodes[1] = &node1;
}

ShellEdge::~ShellEdge()
{
	for (auto &node : nodes)
	{
		node->inclInEdges.remove(this);
	}
	for (auto attedNode : attachedNodes)
	{
		attedNode->isShellEdgeFixed = false;
		attedNode->belongsToShellEdge = nullptr;
	}
}