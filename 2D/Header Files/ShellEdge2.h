#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;
using std::list;

class ShellEdge2 : public unique_ptr_helper<ShellEdge2>
{
public:
	unique_ptr<ShellNode2>* nodes[2];
	vector<unique_ptr<Node2>*> attachedNodes;
	vector<Vector2> attachedNodesStartVectorsToEdge;
	size_t inclInCrysesNum = 0;

	const double Magnitude() const;
	const double SqrMagnitude() const;

	void AttachNodes(const vector<unique_ptr<Node2>*>& free_nodes);
	void ChangeAttachedNode(const size_t& index, const vector<unique_ptr<Node2>*>& free_nodes);
	void SetAttachedNodesStartVectorsToEdge();
	void SetAttachedNodesDistanceFromStartPositionToEdge(const double& units, const double& outOf);

	const bool ContainsAttachedNode(const unique_ptr<Node2>* const& node);
	const bool IsContaining(const ShellNode2& node) const;

	ShellEdge2();
	ShellEdge2(ShellNode2& node0, ShellNode2& node1);
	~ShellEdge2();
};