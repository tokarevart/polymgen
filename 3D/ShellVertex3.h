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
	unique_ptr<Vector3> _position;

public:
	vector<ShellEdge3*> inclInEdges;

	Vector3& GetPosition() const;

	unique_ptr<Vertex3>* FindAttachedVertex(const vector<unique_ptr<Vertex3>*>& free_nodes);

	double& operator[](const int& axisIndex);
	Vector3 operator-(const ShellVertex3& right) const;
	Vector3 operator-(const Vertex3& right) const;

	ShellVertex3();
	ShellVertex3(const double& coor0, const double& coor1, const double& coor2);
	ShellVertex3(const Vector3& position);
	~ShellVertex3();
};