#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Simplex2
{
public:
	Edge* edges[3];
	Crystallite2* inclInCrys;

	bool IsContaining(Edge &edge);
	bool IsContaining(Node2 &node);

	Simplex2();
	Simplex2(Edge edge0, Edge edge1, Edge edge2);
	~Simplex2();
};