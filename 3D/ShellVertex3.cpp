#include "ShellVertex3.h"


Vector3& ShellVertex3::GetPosition() const
{
	return *_position;
}

unique_ptr<Vertex3>* ShellVertex3::FindAttachedVertex(const vector<unique_ptr<Vertex3>*>& freeNodes)
{
	for (auto &vertex : freeNodes)
		if ((*vertex)->belongsToShellVertex == this)
			return vertex;

	return nullptr;
}

double& ShellVertex3::operator[](const int& axisIndex)
{
	return (*_position)[axisIndex];
}

Vector3 ShellVertex3::operator-(const ShellVertex3& right) const
{
	return *_position - *right._position;
}

Vector3 ShellVertex3::operator-(const Vertex3& right) const
{
	return *_position - right.GetPosition();
}

ShellVertex3::ShellVertex3()
{
	_position.reset(new Vector3());
}

ShellVertex3::ShellVertex3(const double& coor0, const double& coor1, const double& coor2)
{
	_position.reset(new Vector3(coor0, coor1, coor2));
}

ShellVertex3::ShellVertex3(const Vector3& position)
{
	_position.reset(new Vector3(position));
}

ShellVertex3::~ShellVertex3() {}