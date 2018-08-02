#include "Node2.h"


void Node2::DestroyIfNoLinks()
{
	if (inclInEdges.empty())
	{
		delete this;
	}
}

Node2::Node2()
{
}

Node2::Node2(double coor0, double coor1)
{
	coors[0] = coor0;
	coors[1] = coor1;
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
}