#include "ShellNode2.h"


const Vector2& ShellNode2::GetPosition()
{
	return *_position;
}

const double& ShellNode2::operator[](const int& axisIndex)
{
	return (*_position)[axisIndex];
}

ShellNode2::ShellNode2()
{
}

ShellNode2::ShellNode2(const ShellNode2& node)
{
	*this = node;
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
	delete _position;
}