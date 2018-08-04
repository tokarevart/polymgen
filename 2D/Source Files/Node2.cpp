#include "Node2.h"


void Node2::SetPosition(const Vector2& newPos)
{
	*_position = newPos;
}

void Node2::SetPosition(const double& coor0, const double& coor1)
{
	(*_position)[0] = coor0;
	(*_position)[1] = coor1;
}

const Vector2& Node2::GetPosition()
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

double& Node2::operator[](const int& axisIndex)
{
	return (*_position)[axisIndex];
}

Node2::Node2()
{
}

Node2::Node2(const Node2& node)
{
	*this = node;
}

Node2::Node2(const double& coor0, const double& coor1)
{
	(*_position)[0] = coor0;
	(*_position)[1] = coor1;
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