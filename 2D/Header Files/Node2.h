#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Node2
{
private:
	Vector2* _position;

public:
	bool isAddedToNodesData = false;
	size_t globalNum;
	ShellEdge2* belongsToShellEdge;
	ShellNode2* belongsToShellNode;
	list<Node2*> neighbors;
	list<Edge2*> inclInEdges;
	list<Simplex2*> inclInSimplexes;
	Polycrystalline2* inclInPolycrys;

	void SetPosition(const Vector2& newPos);
	const const Vector2 GetPosition();
	void DestroyIfNoLinks();

	Node2& operator=(const Node2& node);

	Node2();
	Node2(const Node2& node);
	Node2(double coor0, double coor1);
	Node2(const Vector2& position);
	~Node2();
};