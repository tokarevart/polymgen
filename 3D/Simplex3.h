#pragma once
#include <algorithm>
#include <memory>
#include "ClassDefinitions.h"
#include "ClassInclusions.h"
#include "unique_ptr_helper.h"

using std::unique_ptr;

class Simplex3 : public unique_ptr_helper<Simplex3>
{
public:
	unique_ptr<Vertex3>* vertexes[4];

	double Volume() const;
	double Quality() const;

	Simplex3();
	Simplex3(
		Vertex3 &vert0,
		Vertex3 &vert1,
		Vertex3 &vert2,
		Vertex3 &vert3);
	~Simplex3();
};