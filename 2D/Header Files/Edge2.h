#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Edge2
{
public:
	std::unique_ptr<Node2>* nodes[2];

	std::list<std::unique_ptr<Simplex2>*> inclInSimplexes;

	const double Magnitude() const;
	const double SqrMagnitude() const;

	const bool IsContaining(const Node2& node);

	void DestroyIfNoLinks();

	Edge2();
	Edge2(Node2& node0, Node2& node1);
	~Edge2();
};