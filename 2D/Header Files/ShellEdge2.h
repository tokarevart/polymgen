#pragma once
#include <list>
#include <memory>
#include "AllClassDefinitions.h"
#include "AllClassInclusions.h"

using std::unique_ptr;
using std::list;

#include "unique_ptr_helper.h"
class ShellEdge2 : public unique_ptr_helper<ShellEdge2>
{
private:
	unique_ptr<Vector2> _normal;

	void SetNormal();

public:
	unique_ptr<ShellNode2>* nodes[2];
	list<Crystallite2*> inclInCryses;

	const double Magnitude() const;
	const double SqrMagnitude() const;

	const bool IsContaining(const ShellNode2& node) const;

	ShellEdge2();
	ShellEdge2(ShellNode2& node0, ShellNode2& node1);
	~ShellEdge2();
};