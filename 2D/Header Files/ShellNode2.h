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
	list<Crystallite2*> inclInCryses;

	const Vector2& GetPosition();

	double& operator[](const int& axisIndex);
	Vector2 operator-(const ShellNode2& right);

	ShellNode2();
	ShellNode2(double coor0, double coor1);
	~ShellNode2();
};