#pragma once
#include <list>
#include <memory>
#include "AllClassDefinitions.h"
#include "AllClassInclusions.h"

using std::unique_ptr;
using std::list;

#include "unique_ptr_helper.h"
class Edge2 : public unique_ptr_helper<Edge2>
{
public:
	unique_ptr<Node2>* nodes[2];

	list<unique_ptr<Simplex2>*> inclInSimplexes;

	const double Magnitude() const;
	const double SqrMagnitude() const;

	const bool IsContaining(const Node2& node);

	void DestroyIfNoLinks();

	Edge2();
	Edge2(Node2& node0, Node2& node1);
	~Edge2();
};