#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class ShellEdge2
{
private:
	Vector2* const _normal = new Vector2();

	void SetNormal();

public:
	ShellNode2* nodes[2];
	list<Node2*> attachedNodes;
	list<Crystallite2*> inclInCryses;

	const double Magnitude() const;
	const double SqrMagnitude() const;
	const bool IsContaining(const ShellNode2& node) const;
	const bool IsAttached(const Node2& node) const;

	ShellEdge2();
	ShellEdge2(ShellNode2& node0, ShellNode2& node1);
	~ShellEdge2();
};