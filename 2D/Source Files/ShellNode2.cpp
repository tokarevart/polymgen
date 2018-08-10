#include "ShellNode2.h"


const Vector2& ShellNode2::GetPosition() const
{
	return *_position;
}

void ShellNode2::DestroyIfNoLinks()
{
	if (inclInEdges.empty())
	{
		delete _uniquePtr->release();
	}
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
	return *_position - right.GetPosition();
}

ShellNode2::ShellNode2() : unique_ptr_helper<ShellNode2>(this)
{
	_position.reset(new Vector2());
}

ShellNode2::ShellNode2(const double& coor0, const double& coor1) : unique_ptr_helper<ShellNode2>(this)
{
	_position.reset(new Vector2(coor0, coor1));
}

ShellNode2::ShellNode2(const Vector2& position) : unique_ptr_helper<ShellNode2>(this)
{
	_position.reset(new Vector2(position));
}

ShellNode2::~ShellNode2()
{
	for (auto &edge : inclInEdges)
	{
		if (*edge)
		{
			delete edge->release();
		}
	}
}