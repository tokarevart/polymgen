#pragma once
#include <list>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;
using std::list;

class Vertex3 : public unique_ptr_helper<Vertex3>
{
private:
	unique_ptr<Vector3> _position;

public:
	size_t globalNum;

	ShellFacet3* belongsToShellFacet = nullptr;
	ShellEdge3* belongsToShellEdge = nullptr;
	ShellVertex3* belongsToShellVertex = nullptr;

	//list<unique_ptr<Vertex3>*> neighbors;
	//list<unique_ptr<Edge3>*> inclInEdges;
	//list<unique_ptr<Facet3>*> inclInFacets;

	void SetPosition(
		const Vector3 &newPos);
	void SetPosition(
		const double &coor0, 
		const double &coor1, 
		const double &coor2);

	const Vector3& GetPosition() const;

	//void DestroyIfNoLinks();

	double& operator[](const int &axisIndex);
	Vector3 operator-(const Vertex3 &right) const;
	Vector3 operator-(const ShellVertex3 &right) const;
	Vertex3& operator+=(const Vector3 &right);
	Vertex3& operator-=(const Vector3 &right);

	Vertex3();
	Vertex3(
		const double &coor0, 
		const double &coor1, 
		const double &coor2);
	Vertex3(
		const Vector3 &position);
	~Vertex3();
};