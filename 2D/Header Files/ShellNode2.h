#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class ShellNode2
{
public:
	double coors[2];
	list<ShellEdge*> inclInEdges;

	ShellNode2();
	ShellNode2(double coor0, double coor1);
	~ShellNode2();
};