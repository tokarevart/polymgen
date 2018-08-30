#include "Polycrystalline3.h"
#include <algorithm>
#include <iostream>


//ShellEdge3* Polycrystalline3::FindShellEdge(ShellVertex3& vertex0, ShellVertex3& vertex1)
//{
//	return nullptr;
//}

unique_ptr<Edge3>* Polycrystalline3::FindEdge(Vertex3& vertex0, Vertex3& vertex1)
{
	if (std::find(vertex0.neighbors.begin(), vertex0.neighbors.end(), vertex1.GetPtrToUniquePtr()) == vertex0.neighbors.end())
	{
		return nullptr;
	}

	for (list<unique_ptr<Edge3>*>::const_iterator edges_iter = vertex0.inclInEdges.begin(), end = vertex0.inclInEdges.end();
		edges_iter != end; 
		edges_iter++)
	{
		if ((**edges_iter)->vertexes[0]->get() == &vertex0 && 
			(**edges_iter)->vertexes[1]->get() == &vertex1 ||
			(**edges_iter)->vertexes[0]->get() == &vertex1 && 
			(**edges_iter)->vertexes[1]->get() == &vertex0)
		{
			return *edges_iter;
		}
	}

	throw std::exception("Error in function Polycrystalline3::FindEdge");
	return nullptr;
}

template <class T>
void ErasePtrsToNullptr(vector<unique_ptr<T>*>& vec)
{
	size_t real_objs_num = std::count_if(vec.begin(), vec.end(), [](unique_ptr<T>*& ptr) { return *ptr ? true : false; });
	vector<unique_ptr<T>*> buf_vec(real_objs_num);

	size_t firts_thr_nums = real_objs_num / 2;
	size_t second_thr_nums = real_objs_num - firts_thr_nums;
	//#pragma omp parallel num_threads(2)
	//{
	//#pragma omp single
	//{
	size_t index1 = 0;
	for (size_t i = 0; index1 < firts_thr_nums; i++)
	{
		if (*vec[i])
		{
			buf_vec[index1] = vec[i];
			index1++;
		}
	}
	//}
	//#pragma omp single
	//{
	size_t index2 = 0;
	for (size_t i = vec.size() - 1; index2 < second_thr_nums; i--)
	{
		if (*vec[i])
		{
			buf_vec[real_objs_num - 1 - index2] = vec[i];
			index2++;
		}
	}
	//}
	//}

	vec = std::move(buf_vec);
}

void Polycrystalline3::ErasePtrsToNullptrFromVectors()
{
	//#pragma omp parallel num_threads(3)
	//{
	//#pragma omp single
	ErasePtrsToNullptr(_freeVertexes);
	//#pragma omp single
	ErasePtrsToNullptr(_freeEdges);
	//#pragma omp single
	ErasePtrsToNullptr(_freeFacets);
	//}
}

void Polycrystalline3::TriangulateShell()
{
	size_t edges_num = _freeEdges.size();
	for (size_t i = 0; i < edges_num; i++)
	{
		if (*_freeEdges[i] &&
			(*_freeEdges[i])->SqrMagnitude() > 2.25 * preferredLength * preferredLength)
		{
			(*_freeEdges[i])->MakeTwoInstead(_freeFacets, _freeEdges, _freeVertexes, _shellEdges);
			edges_num += 2;
		}
	}
}

Vector3 Polycrystalline3::ShiftToLavEdges(const Vertex3& vertex)
{
	Vector3 shift;
	const double b = 0.15;
	for (auto &neighbor : vertex.neighbors)
	{
		if (*neighbor)
		{
			Vector3 vec = **neighbor - vertex;
			double buf = b * (1.0 - preferredLength / vec.Magnitude());
			if ((*neighbor)->belongsToShellVertex)
			{
				shift += vec * 2.0 * buf;
			}
			else if ((*neighbor)->belongsToShellEdge)
			{
				double k = 2.0 - abs(Vector3::Cos(vec, *(*neighbor)->belongsToShellEdge->vertexes[0] - *(*neighbor)->belongsToShellEdge->vertexes[1]));
				shift += vec * k * buf;
			}
			else
			{
				shift += vec * buf;
			}
		}
	}

	if (shift.SqrMagnitude() > preferredLength * preferredLength * 0.0025)
	{
		shift = shift.Normalize() * preferredLength * 0.05;
	}

	return shift;
}

void Polycrystalline3::DistributeVertexesOnShellEvenly(size_t iterations_num)
{
	size_t verts_num = _freeVertexes.size();
	Vector3* shifts = new Vector3[verts_num];

	if (iterations_num == 0)
	{
		std::cout << "Enter the iterations number: ";
		std::cin >> iterations_num;
	}
	for (int i = 0; i < iterations_num; i++)
	{
		//#pragma omp parallel firstprivate(verts_num)
		//{
		//#pragma omp for
		for (size_t j = 0; j < verts_num; j++)
		{
			if (*_freeVertexes[j])
			{
				shifts[j] = ShiftToLavEdges(**_freeVertexes[j]);
			}
		}
		//#pragma omp for
		for (size_t j = 0; j < verts_num; j++)
		{
			if (*_freeVertexes[j])
			{
				**_freeVertexes[j] += shifts[j];

				//if ((*_freeVertexes[j])->belongsToShellEdge)
				//{
				//	(*_freeVertexes[j])->SetPosition(
				//		Vector3::Project(
				//		(*_freeVertexes[j])->GetPosition(),
				//			(*_freeVertexes[j])->belongsToShellEdge->vertexes[0]->GetPosition(),
				//			(*_freeVertexes[j])->belongsToShellEdge->vertexes[1]->GetPosition()));
				//}
				//else if ((*_freeVertexes[j])->belongsToShellFacet)
				//{
				//	Vector3 plane_point, a, b;
				//	plane_point = (*_freeVertexes[j])->belongsToShellFacet->edges[0]->vertexes[0]->GetPosition();
				//	a = (*_freeVertexes[j])->belongsToShellFacet->edges[0]->vertexes[1]->GetPosition() - plane_point;
				//	if ((*_freeVertexes[j])->belongsToShellFacet->edges[0]->vertexes[0] == (*_freeVertexes[j])->belongsToShellFacet->edges[1]->vertexes[0])
				//	{
				//		b = (*_freeVertexes[j])->belongsToShellFacet->edges[1]->vertexes[1]->GetPosition() - plane_point;
				//	}
				//	else
				//	{
				//		b = (*_freeVertexes[j])->belongsToShellFacet->edges[1]->vertexes[0]->GetPosition() - plane_point;
				//	}

				//	(*_freeVertexes[j])->SetPosition(
				//		plane_point + ((*_freeVertexes[j])->GetPosition() - plane_point).Project(a, b));
				//}
			}
		}
		//}
	}

	delete[] shifts;
}

void Polycrystalline3::DelaunayPostprocessing()
{
	unique_ptr<Vertex3>* around_nodes[2];
	size_t edges_num = _freeEdges.size();
	for (size_t i = 0; i < edges_num; i++) // Can be partialy parallelized
	{
		if (!*_freeEdges[i] ||
			((*(*_freeEdges[i])->vertexes[0])->belongsToShellEdge == (*(*_freeEdges[i])->vertexes[1])->belongsToShellEdge && 
			 (*(*_freeEdges[i])->vertexes[1])->belongsToShellEdge) ||
			((*(*_freeEdges[i])->vertexes[0])->belongsToShellEdge && 
			 (*(*_freeEdges[i])->vertexes[1])->belongsToShellVertex && 
			 (*(*_freeEdges[i])->vertexes[0])->belongsToShellEdge->IsContaining(*(*(*_freeEdges[i])->vertexes[1])->belongsToShellVertex)) ||
			((*(*_freeEdges[i])->vertexes[1])->belongsToShellEdge &&
			 (*(*_freeEdges[i])->vertexes[0])->belongsToShellVertex &&
			 ((*(*_freeEdges[i])->vertexes[1])->belongsToShellEdge)->IsContaining(*(*(*_freeEdges[i])->vertexes[0])->belongsToShellVertex)))
		{
			continue;
		}
		
		around_nodes[0] = (*(*_freeEdges[i])->inclInFacets.front())->FindVertexNotIncludedInEdge(**_freeEdges[i]);
		around_nodes[1] = (*(*_freeEdges[i])->inclInFacets.back())->FindVertexNotIncludedInEdge(**_freeEdges[i]);
		if ((*around_nodes[0])->belongsToShellEdge && (*around_nodes[1])->belongsToShellVertex ||
			(*around_nodes[1])->belongsToShellEdge && (*around_nodes[0])->belongsToShellVertex)
		{
			continue;
		}

		(*_freeEdges[i])->FlipIfNeeded(_freeEdges, _freeFacets);
	}
}

void Polycrystalline3::InputData()
{
	preferredLength = 0.1;

	_shellVertexes.insert(_shellVertexes.end(), 
	{
		new ShellVertex3(1.0, 0.0, 0.0),
		new ShellVertex3(0.5, 0.0, 1.0),
		new ShellVertex3(0.0, 0.0, 0.0),
		new ShellVertex3(0.5, 1.0, 0.0)
	});
	_shellEdges.insert(_shellEdges.end(),
	{
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[1]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[2]),
		new ShellEdge3(*_shellVertexes[2], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[2])
	});
	_shellFacets.insert(_shellFacets.end(),
	{
		new ShellFacet3(*_shellEdges[0], *_shellEdges[1], *_shellEdges[5]),
		new ShellFacet3(*_shellEdges[1], *_shellEdges[2], *_shellEdges[3]),
		new ShellFacet3(*_shellEdges[0], *_shellEdges[3], *_shellEdges[4]),
		new ShellFacet3(*_shellEdges[2], *_shellEdges[4], *_shellEdges[5])
	});

	crystallites.push_back(new Crystallite3);
	crystallites.front()->shellEdges = _shellEdges;
	crystallites.front()->shellFacets = _shellFacets;


	_freeVertexes.insert(_freeVertexes.end(),
	{
		(new Vertex3(1.0, 0.0, 0.0))->GetPtrToUniquePtr(),
		(new Vertex3(0.5, 0.0, 1.0))->GetPtrToUniquePtr(),
		(new Vertex3(0.0, 0.0, 0.0))->GetPtrToUniquePtr(),
		(new Vertex3(0.5, 1.0, 0.0))->GetPtrToUniquePtr()
	});
	for (size_t i = 0, shell_vertexes_num = _shellVertexes.size(); i < shell_vertexes_num; i++)
	{
		(*_freeVertexes[i])->belongsToShellVertex = _shellVertexes[i];
	}
	_freeEdges.insert(_freeEdges.end(),
	{
		(new Edge3(**_freeVertexes[0], **_freeVertexes[1]))->GetPtrToUniquePtr(),
		(new Edge3(**_freeVertexes[1], **_freeVertexes[2]))->GetPtrToUniquePtr(),
		(new Edge3(**_freeVertexes[2], **_freeVertexes[3]))->GetPtrToUniquePtr(),
		(new Edge3(**_freeVertexes[1], **_freeVertexes[3]))->GetPtrToUniquePtr(),
		(new Edge3(**_freeVertexes[0], **_freeVertexes[3]))->GetPtrToUniquePtr(),
		(new Edge3(**_freeVertexes[0], **_freeVertexes[2]))->GetPtrToUniquePtr()
	});
	_freeFacets.insert(_freeFacets.end(),
	{
		(new Facet3(**_freeEdges[0], **_freeEdges[1], **_freeEdges[5]))->GetPtrToUniquePtr(),
		(new Facet3(**_freeEdges[1], **_freeEdges[2], **_freeEdges[3]))->GetPtrToUniquePtr(),
		(new Facet3(**_freeEdges[0], **_freeEdges[3], **_freeEdges[4]))->GetPtrToUniquePtr(),
		(new Facet3(**_freeEdges[2], **_freeEdges[4], **_freeEdges[5]))->GetPtrToUniquePtr()
	});
}

void Polycrystalline3::OutputData(string filename) const
{
	ofstream file(filename);
	
	for (size_t i = 0, verts_num = _freeVertexes.size(); i < verts_num; i++)
	{
		file << "v " << (**_freeVertexes[i])[0] << ' ' << (**_freeVertexes[i])[1] << ' ' << (**_freeVertexes[i])[2] << '\n';
		(*_freeVertexes[i])->globalNum = i + 1;
	}
	for (auto &facet : _freeFacets)
	{
		if (!*facet)
		{
			continue;
		}

		vector<size_t> gl_nums;
		for (auto &edge : (*facet)->edges)
		{
			for (auto &vert : (*edge)->vertexes)
			{
				if (std::find(gl_nums.begin(), gl_nums.end(), (*vert)->globalNum) == gl_nums.end())
				{
					gl_nums.push_back((*vert)->globalNum);
				}
			}
		}
		file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
	}

	file.close();
}

Polycrystalline3::Polycrystalline3() {}

Polycrystalline3::~Polycrystalline3()
{
	// Delete simplexes

	for (auto &facet : _shellFacets)
	{
		delete facet;
	}
	for (auto &edge : _shellEdges)
	{
		delete edge;
	}
	for (auto &vert : _shellVertexes)
	{
		delete vert;
	}

	for (auto &ptr : _freeFacets)
	{
		delete ptr;
	}
	for (auto &ptr : _freeEdges)
	{
		delete ptr;
	}
	for (auto &ptr : _freeVertexes)
	{
		delete ptr;
	}
}