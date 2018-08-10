#include "ShellEdge2.h"

const double MINUS_PI_DIV_2 = -3.141592653589793 * 0.5;

void ShellEdge2::SetNormal()
{
	(*_normal)[0] = (*nodes[1]->get())[0] - (*nodes[0]->get())[0];
	(*_normal)[1] = (*nodes[1]->get())[1] - (*nodes[0]->get())[1];

	_normal->Rotate(MINUS_PI_DIV_2, Radian).Normalize();
}

const double ShellEdge2::Magnitude() const
{
	return sqrt(SqrMagnitude());
}

const double ShellEdge2::SqrMagnitude() const
{
	Vector2 buf = **nodes[1] - **nodes[0];
	return Vector2::DotProduct(buf, buf);
}

const bool ShellEdge2::IsContaining(const ShellNode2& node) const
{
	if (nodes[0]->get() == &node ||
		nodes[1]->get() == &node)
	{
		return true;
	}

	return false;
}

ShellEdge2::ShellEdge2() : unique_ptr_helper<ShellEdge2>(this)
{
	nodes[0] = nullptr;
	nodes[1] = nullptr;
	_normal.reset(new Vector2());
}

ShellEdge2::ShellEdge2(ShellNode2& node0, ShellNode2& node1) : unique_ptr_helper<ShellEdge2>(this)
{
	_normal.reset(new Vector2());

	nodes[0] = node0.GetPtrToUniquePtr();
	nodes[1] = node1.GetPtrToUniquePtr();

	(*nodes[0])->inclInEdges.push_back(GetPtrToUniquePtr());
	(*nodes[1])->inclInEdges.push_back(GetPtrToUniquePtr());

	SetNormal();
}

ShellEdge2::~ShellEdge2()
{
	for (auto &node : nodes)
	{
		if (*node)
		{
			(*node)->inclInEdges.remove(GetPtrToUniquePtr());
			(*node)->DestroyIfNoLinks();
		}
	}
}