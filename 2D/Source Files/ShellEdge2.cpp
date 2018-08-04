#include "ShellEdge2.h"


bool ShellEdge2::Contains(ShellNode2 &node)
{
	if (nodes[0] == &node ||
		nodes[1] == &node)
	{
		return true;
	}

	return false;
}

bool ShellEdge2::IsAttached(Node2 &node)
{
	return find(attachedNodes.begin(), attachedNodes.end(), &node) != attachedNodes.end();
}

ShellEdge2::ShellEdge2()
{
}

ShellEdge2::ShellEdge2(ShellNode2 &node0, ShellNode2 &node1)
{
	nodes[0] = &node0;
	nodes[1] = &node1;

	// Calculate normal vector.
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
}