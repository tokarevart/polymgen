#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;
using std::list;

class ShellNode2 : public unique_ptr_helper<ShellNode2>
{
private:
	unique_ptr<Vector2> _position;

public:
	list<unique_ptr<ShellEdge2>*> inclInEdges;
	list<Crystallite2*> inclInCryses;

	const Vector2& GetPosition() const;

	void DestroyIfNoLinks();

	double& operator[](const int& axisIndex);
	Vector2 operator-(const ShellNode2& right) const;
	Vector2 operator-(const Node2& right) const;

	ShellNode2();
	ShellNode2(const double& coor0, const double& coor1);
	ShellNode2(const Vector2& position);
	~ShellNode2();
};