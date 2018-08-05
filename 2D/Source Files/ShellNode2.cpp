#include "ShellNode2.h"


const Vector2& ShellNode2::GetPosition()
{
	return *_position;
}

double& ShellNode2::operator[](const int& axisIndex)
{
	return (*_position)[axisIndex];
}

Vector2 ShellNode2::operator-(const ShellNode2& right)
{
	return *_position - *right._position;
}

ShellNode2::ShellNode2()
{
}

ShellNode2::ShellNode2(const double& coor0, const double& coor1)
{
	(*_position)[0] = coor0;
	(*_position)[1] = coor1;
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