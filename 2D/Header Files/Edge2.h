#pragma once
#include <list>
#include <vector>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;
using std::list;
using std::vector;

class Edge2 : public unique_ptr_helper<Edge2>
{
public:
	unique_ptr<Node2>* nodes[2];

	list<unique_ptr<Simplex2>*> inclInSimplexes;

	const double Magnitude() const;
	const double SqrMagnitude() const;

	void MakeTwoInstead(list<unique_ptr<Simplex2>*>& freeSimplexes, vector<unique_ptr<Edge2>*>& freeEdges, vector<unique_ptr<Node2>*>& freeNodes);

	const bool IsContaining(const Node2& node);
	const bool BelongsToShell();
	const bool IntersectsWith(const Edge2& edge);
	const bool IntersectsWith(const ShellEdge2& edge);
	const bool IntersectsWith(const vector<unique_ptr<ShellEdge2>*>& shellEdges);

	void DestroyIfNoLinks();

	Edge2();
	Edge2(Node2& node0, Node2& node1);
	~Edge2();
};