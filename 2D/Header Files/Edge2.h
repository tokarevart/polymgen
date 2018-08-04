#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Edge2
{
public:
	Node2* nodes[2];
	list<Simplex2*> inclInSimplexes;

	const bool IsContaining(const Node2& node);
	void DestroyIfNoLinks();

	Edge2();
	Edge2(Node2& node0, Node2& node1);
	~Edge2();
};