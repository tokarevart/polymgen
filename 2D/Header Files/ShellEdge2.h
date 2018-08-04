#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class ShellEdge2
{
public:
	ShellNode2* nodes[2];
	Vector2* normal;
	list<Node2*> attachedNodes;

	bool Contains(ShellNode2 &node);
	bool IsAttached(Node2 &node);

	ShellEdge2();
	ShellEdge2(ShellNode2 &node0, ShellNode2 &node1);
	~ShellEdge2();
};