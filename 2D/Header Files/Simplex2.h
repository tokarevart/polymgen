#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Simplex2
{
public:
	Edge2* edges[3];
	Crystallite2* inclInCrys;

	bool IsContaining(Edge2 &edge);
	bool IsContaining(Node2 &node);

	Simplex2();
	Simplex2(Edge2 &edge0, Edge2 &edge1, Edge2 &edge2);
	Simplex2(Node2 &node0, Node2 &node1, Node2 &node2);
	~Simplex2();
};