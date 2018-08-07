#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Node2
{
private:
	Vector2* _position;

public:
	bool isAddedToNodesData = false;
	size_t globalNum;

	ShellEdge2* belongsToShellEdge;
	ShellNode2* belongsToShellNode;

	std::list<Node2*> neighbors;
	std::list<Edge2*> inclInEdges;
	std::list<Simplex2*> inclInSimplexes;
	std::list<Crystallite2*> inclInCryses;

	void SetPosition(const Vector2& newPos);
	void SetPosition(const double& coor0, const double& coor1);

	const Vector2& GetPosition() const;

	void DestroyIfNoLinks();

	double& operator[](const int& axisIndex);
	Vector2 operator-(const Node2& right) const;
	Node2& operator+=(const Vector2& right);
	Node2& operator-=(const Vector2& right);

	Node2();
	Node2(const double& coor0, const double& coor1);
	Node2(const Vector2& position);
	~Node2();
};