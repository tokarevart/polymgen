#include "Node2.h"


void Node2::SetPosition(const Vector2& newPos)
{
	*_position = newPos;
}

const Vector2 Node2::GetPosition()
{
	return *_position;
}

void Node2::DestroyIfNoLinks()
{
	if (inclInEdges.empty())
	{
		delete this;
	}
}

Node2::Node2()
{
	_position = new Vector2();
}

Node2::Node2(const Node2& node)
{
	*this = node;
}

Node2::Node2(double coor0, double coor1)
{
	_position = new Vector2(coor0, coor1);
}

Node2::Node2(const Vector2& position)
{
	*_position = position;
}

Node2::~Node2()
{
	for (auto &edge : inclInEdges)
	{
		if (edge)
		{
			delete edge;
		}
	}
	delete _position;
}