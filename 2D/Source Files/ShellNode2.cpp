#include "ShellNode2.h"


const Vector2& ShellNode2::GetPosition() const
{
	return *_position;
}

double& ShellNode2::operator[](const int& axisIndex)
{
	return (*_position)[axisIndex];
}

Vector2 ShellNode2::operator-(const ShellNode2& right) const
{
	return *_position - *right._position;
}

Vector2 ShellNode2::operator-(const Node2& right) const
{
	return *_position - *right._position;
}

ShellNode2::ShellNode2() 
{
	_position = new Vector2();
}

ShellNode2::ShellNode2(const double& coor0, const double& coor1)
{
	_position = new Vector2(coor0, coor1);
}

ShellNode2::~ShellNode2()
{
	for (auto &edge : inclInEdges)
	{
		if (edge)
		{
			delete edge;
		}
	}
	if (_position)
	{
		delete _position;
	}
}