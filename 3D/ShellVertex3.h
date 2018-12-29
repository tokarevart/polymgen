#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"

using std::unique_ptr;
using std::list;

class ShellVertex3
{
private:
	unique_ptr<Vec3> _position;

public:
	//vector<ShellEdge3*> inclInEdges;

	Vec3& getPosition() const;

	unique_ptr<Vertex3>* FindAttachedVertex(
		const vector<unique_ptr<Vertex3>*> &freeNodes);

	double& operator[](const int &axisIndex);
	Vec3 operator-(const ShellVertex3 &right) const;
	Vec3 operator-(const Vertex3 &right) const;

	ShellVertex3();
	ShellVertex3(
		const double &coor0, 
		const double &coor1, 
		const double &coor2);
	ShellVertex3(
		const Vec3 &position);
	~ShellVertex3();
};