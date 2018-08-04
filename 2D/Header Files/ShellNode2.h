#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class ShellNode2
{
private:
	Vector2* _position;

public:
	Node2* attachedNode;
	list<ShellEdge2*> inclInEdges;

	const const Vector2 GetPosition();

	Node2& operator=(const Node2& node);

	ShellNode2();
	ShellNode2(const ShellNode2& node);
	ShellNode2(double coor0, double coor1);
	~ShellNode2();
};