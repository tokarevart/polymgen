#pragma once
#include <list>
#include <memory>
#include "AllClassDefinitions.h"
#include "AllClassInclusions.h"

using std::unique_ptr;
using std::list;

#include "unique_ptr_helper.h"
class ShellNode2 : public unique_ptr_helper<ShellNode2>
{
private:
	unique_ptr<Vector2> _position;

public:
	list<unique_ptr<ShellEdge2>*> inclInEdges;
	list<Crystallite2*> inclInCryses;

	const Vector2& GetPosition() const;

	double& operator[](const int& axisIndex);
	Vector2 operator-(const ShellNode2& right) const;
	Vector2 operator-(const Node2& right) const;

	ShellNode2();
	ShellNode2(const double& coor0, const double& coor1);
	ShellNode2(const Vector2& position);
	~ShellNode2();
};