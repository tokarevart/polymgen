#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Edge
{
public:
	bool isShellPart = false;
	Node2* nodes[2];
	list<Simplex2*> inclInSimplexes;

	bool IsContaining(Node2 &node);
	void DestroyIfNoLinks();

	Edge();
	Edge(Node2 &node0, Node2 &node1);
	~Edge();
};