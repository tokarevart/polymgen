#pragma once
#include <list>
#include "AllClassInclusions.h"
#include "AllClassDefinitions.h"

using namespace std;


class Node2
{
public:
	bool isShellEdgeFixed = false;
	bool isShellNodeFixed = false;
	bool isAddedToNodesData = false;
	double coors[2];
	size_t globalNum;
	ShellEdge* belongsToShellEdge;
	ShellNode2* belongsToShellNode;
	list<Node2*> neighbors;
	list<Edge*> inclInEdges;
	list<Simplex2*> inclInSimplexes;
	Polycrystalline* inclInPolycrys;

	void DestroyIfNoLinks();

	Node2();
	Node2(double coor0, double coor1);
	~Node2();
};