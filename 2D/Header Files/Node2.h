#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"


class Node2
{
private:
	std::unique_ptr<Vector2> _position;

public:
	bool isAddedToNodesData = false;
	size_t globalNum;

	std::unique_ptr<ShellEdge2>* belongsToShellEdge = nullptr;
	std::unique_ptr<ShellNode2>* belongsToShellNode = nullptr;

	std::list<std::unique_ptr<Node2>*> neighbors;
	std::list<std::unique_ptr<Edge2>*> inclInEdges;
	std::list<std::unique_ptr<Simplex2>*> inclInSimplexes;

	void SetPosition(const Vector2& newPos);
	void SetPosition(const double& coor0, const double& coor1);

	const Vector2& GetPosition() const;

	void DestroyIfNoLinks();

	double& operator[](const int& axisIndex);
	Vector2 operator-(const Node2& right) const;
	Vector2 operator-(const ShellNode2& right) const;
	Node2& operator+=(const Vector2& right);
	Node2& operator-=(const Vector2& right);

	Node2();
	Node2(const double& coor0, const double& coor1);
	Node2(const Vector2& position);
	~Node2();

	friend ShellNode2;
};