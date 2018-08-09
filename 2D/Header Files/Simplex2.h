#pragma once
#include <list>
#include <memory>
#include "AllClassDefinitions.h"
#include "AllClassInclusions.h"

using std::unique_ptr;

#include "unique_ptr_helper.h"
class Simplex2 : public unique_ptr_helper<Simplex2>
{
public:
	unique_ptr<Edge2>* edges[3];

	const bool IsContaining(const Edge2& edge) const;
	const bool IsContaining(const Node2& node) const;

	Simplex2();
	Simplex2(Edge2& edge0, Edge2& edge1, Edge2& edge2);
	Simplex2(Node2& node0, Node2& node1, Node2& node2);
	~Simplex2();
};