#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class ShellNode2
{
private:
	Vector2* _position;

public:
	Node2* attachedNode;
	std::list<ShellEdge2*> inclInEdges;
	std::list<Crystallite2*> inclInCryses;

	const Vector2& GetPosition() const;

	double& operator[](const int& axisIndex);
	Vector2 operator-(const ShellNode2& right) const;

	ShellNode2();
	ShellNode2(const double& coor0, const double& coor1);
	~ShellNode2();
};