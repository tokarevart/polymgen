#include "ShellNode2.h"


const Vector2 ShellNode2::GetPosition()
{
	return *_position;
}

ShellNode2::ShellNode2()
{
	_position = new Vector2();
}

ShellNode2::ShellNode2(const ShellNode2& node)
{
	*this = node;
}

ShellNode2::ShellNode2(double coor0, double coor1)
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
	delete _position;
}