#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class ShellNode2
{
private:
	Vector2* const _position = new Vector2();

public:
	Node2* attachedNode;
	list<ShellEdge2*> inclInEdges;

	const Vector2& GetPosition();

	ShellNode2& operator=(const ShellNode2 &node);
	const double& operator[](const int& axisIndex);

	ShellNode2();
	ShellNode2(const ShellNode2 &node);
	ShellNode2(const double &coor0, const double &coor1);
	~ShellNode2();
};