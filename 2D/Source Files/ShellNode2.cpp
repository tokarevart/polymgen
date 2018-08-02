#include "ShellNode2.h"


ShellNode2::ShellNode2()
{
}

ShellNode2::ShellNode2(double coor0, double coor1)
{
	coors[0] = coor0;
	coors[1] = coor1;
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
}