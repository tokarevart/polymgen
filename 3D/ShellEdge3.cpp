#include "ShellEdge3.h"


const double ShellEdge3::Magnitude() const
{
	return sqrt(SqrMagnitude());
}

const double ShellEdge3::SqrMagnitude() const
{
	Vector3 buf = *vertexes[1] - *vertexes[0];
	return Vector3::DotProduct(buf, buf);
}

const bool ShellEdge3::IsContaining(const ShellVertex3& vertex) const
{
	if (vertexes[0] == &vertex ||
		vertexes[1] == &vertex)
	{
		return true;
	}

	return false;
}

ShellEdge3::ShellEdge3()
{
	vertexes[0] = nullptr;
	vertexes[1] = nullptr;
}

ShellEdge3::ShellEdge3(ShellVertex3& vertex0, ShellVertex3& vertex1)
{
	vertexes[0] = &vertex0;
	vertexes[1] = &vertex1;

	vertexes[0]->inclInEdges.push_back(this);
	vertexes[1]->inclInEdges.push_back(this);
}

ShellEdge3::~ShellEdge3()
{
}