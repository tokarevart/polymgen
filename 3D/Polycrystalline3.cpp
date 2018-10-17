#include "Polycrystalline3.h"
#include <algorithm>
#include <iostream>


void Polycrystalline3::SetLinksWithShell()
{
	double sqr_sufficient_dist = _preferredLength * _preferredLength * 1e-6;
	size_t shell_verts_num = _shellVertexes.size();
	size_t shell_edges_num = _shellEdges.size();
	size_t shell_facets_num = _shellFacets.size();
	//#pragma omp parallel firstprivate(sqr_sufficient_dist, shell_verts_num, shell_edges_num, shell_facets_num)
	//{
	//#pragma omp for
	for (size_t i = 0; i < shell_verts_num; i++)
	{
		for (auto &vert : _startFrontVertexes)
		{
			if ((*_shellVertexes[i] - **vert).SqrMagnitude() < sqr_sufficient_dist)
			{
				(*vert)->belongsToShellVertex = _shellVertexes[i];
				break;
			}
		}
	}
	//#pragma omp for
	for (size_t i = 0; i < shell_edges_num; i++)
	{
		Vector3 proj_buf;
		for (auto &vert : _startFrontVertexes)
		{
			if ((*vert)->belongsToShellVertex)
				continue;

			if (Vector3::Project(proj_buf, (*vert)->GetPosition(), _shellEdges[i]->vertexes[0]->GetPosition(), _shellEdges[i]->vertexes[1]->GetPosition()) &&
				(proj_buf - (*vert)->GetPosition()).SqrMagnitude() < sqr_sufficient_dist)
				(*vert)->belongsToShellEdge = _shellEdges[i];
		}
	}
	//#pragma omp for
	for (size_t i = 0; i < shell_facets_num; i++)
	{
		Vector3 proj_buf;
		for (auto &vert : _startFrontVertexes)
		{
			if ((*vert)->belongsToShellVertex ||
				(*vert)->belongsToShellEdge)
				continue;

			proj_buf = _shellFacets[i]->edges[0]->vertexes[0]->GetPosition() + (**vert - *_shellFacets[i]->edges[0]->vertexes[0]).Project(*_shellFacets[i]->edges[0]->vertexes[1] - *_shellFacets[i]->edges[0]->vertexes[0], *_shellFacets[i]->edges[1]->vertexes[1] - *_shellFacets[i]->edges[1]->vertexes[0]);
			if ((proj_buf - (*vert)->GetPosition()).SqrMagnitude() < sqr_sufficient_dist)
				(*vert)->belongsToShellFacet = _shellFacets[i];
		}
	}
	//}
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
	ErasePtrsToNullptr(_startFrontVertexes);
	//#pragma omp single
	ErasePtrsToNullptr(_startFrontEdges);
	//#pragma omp single
	ErasePtrsToNullptr(_startFrontFacets);
	//}
}

void Polycrystalline3::TriangulateShell()
{
	size_t edges_num = _startFrontEdges.size();
	for (size_t i = 0; i < edges_num; i++)
	{
		if (*_startFrontEdges[i] &&
			(*_startFrontEdges[i])->SqrMagnitude() > 2.25 * _preferredLength * _preferredLength)
		{
			(*_startFrontEdges[i])->MakeTwoInstead(_startFrontFacets, _startFrontEdges, _startFrontVertexes);//, _shellEdges, _shellFacets);
			edges_num += 2;
		}
	}
}

//Vector3 Polycrystalline3::ShiftToLavEdges(const Vertex3& vertex)
//{
//	Vector3 shift;
//	const double b = 0.15;
//	for (auto &neighbor : vertex.neighbors)
//	{
//		if (*neighbor)
//		{
//			Vector3 vec = **neighbor - vertex;
//			double buf = b * (1.0 - _preferredLength / vec.Magnitude());
//			if ((*neighbor)->belongsToShellVertex)
//			{
//				shift += vec * 2.0 * buf;
//			}
//			else if ((*neighbor)->belongsToShellEdge)
//			{
//				double k = 2.0 - abs(Vector3::Cos(vec, *(*neighbor)->belongsToShellEdge->vertexes[0] - *(*neighbor)->belongsToShellEdge->vertexes[1]));
//				shift += vec * k * buf;
//			}
//			else
//			{
//				shift += vec * buf;
//			}
//		}
//	}
//
//	if (shift.SqrMagnitude() > _preferredLength * _preferredLength * 0.0025)
//	{
//		shift = shift.Normalize() * _preferredLength * 0.05;
//	}
//
//	return shift;
//}

//void Polycrystalline3::DistributeVertexesOnShellEvenly(size_t iterations_num)
//{
//	size_t verts_num = _startFrontVertexes.size();
//	Vector3* shifts = new Vector3[verts_num];
//
//	if (iterations_num == 0)
//	{
//		std::cout << "Enter the iterations number: ";
//		std::cin >> iterations_num;
//	}
//	for (int i = 0; i < iterations_num; i++)
//	{
//		//#pragma omp parallel firstprivate(verts_num)
//		//{
//		//#pragma omp for
//		for (size_t j = 0; j < verts_num; j++)
//		{
//			if (*_startFrontVertexes[j])
//			{
//				shifts[j] = ShiftToLavEdges(**_startFrontVertexes[j]);
//			}
//		}
//		//#pragma omp for
//		for (size_t j = 0; j < verts_num; j++)
//		{
//			if (*_startFrontVertexes[j])
//			{
//				**_startFrontVertexes[j] += shifts[j];
//
//				//if ((*_startFrontVertexes[j])->belongsToShellEdge)
//				//{
//				//	(*_startFrontVertexes[j])->SetPosition(
//				//		Vector3::Project(
//				//		(*_startFrontVertexes[j])->GetPosition(),
//				//			(*_startFrontVertexes[j])->belongsToShellEdge->vertexes[0]->GetPosition(),
//				//			(*_startFrontVertexes[j])->belongsToShellEdge->vertexes[1]->GetPosition()));
//				//}
//				//else if ((*_startFrontVertexes[j])->belongsToShellFacet)
//				//{
//				//	Vector3 plane_point, a, b;
//				//	plane_point = (*_startFrontVertexes[j])->belongsToShellFacet->edges[0]->vertexes[0]->GetPosition();
//				//	a = (*_startFrontVertexes[j])->belongsToShellFacet->edges[0]->vertexes[1]->GetPosition() - plane_point;
//				//	if ((*_startFrontVertexes[j])->belongsToShellFacet->edges[0]->vertexes[0] == (*_startFrontVertexes[j])->belongsToShellFacet->edges[1]->vertexes[0])
//				//	{
//				//		b = (*_startFrontVertexes[j])->belongsToShellFacet->edges[1]->vertexes[1]->GetPosition() - plane_point;
//				//	}
//				//	else
//				//	{
//				//		b = (*_startFrontVertexes[j])->belongsToShellFacet->edges[1]->vertexes[0]->GetPosition() - plane_point;
//				//	}
//				//
//				//	(*_startFrontVertexes[j])->SetPosition(
//				//		plane_point + ((*_startFrontVertexes[j])->GetPosition() - plane_point).Project(a, b));
//				//}
//			}
//		}
//		//}
//	}
//
//	delete[] shifts;
//}

void Polycrystalline3::DelaunayPostprocessing()
{
	unique_ptr<Vertex3>* around_nodes[2];
	unique_ptr<Facet3>* around_facets[2];
	size_t edges_num = _startFrontEdges.size();
	for (size_t i = 0; i < edges_num; i++)
	{
		if (!*_startFrontEdges[i] ||
			((*(*_startFrontEdges[i])->vertexes[0])->belongsToShellEdge == (*(*_startFrontEdges[i])->vertexes[1])->belongsToShellEdge && 
			 (*(*_startFrontEdges[i])->vertexes[1])->belongsToShellEdge) ||
			((*(*_startFrontEdges[i])->vertexes[0])->belongsToShellEdge && 
			 (*(*_startFrontEdges[i])->vertexes[1])->belongsToShellVertex && 
			 (*(*_startFrontEdges[i])->vertexes[0])->belongsToShellEdge->IsContaining(*(*(*_startFrontEdges[i])->vertexes[1])->belongsToShellVertex)) ||
			((*(*_startFrontEdges[i])->vertexes[1])->belongsToShellEdge &&
			 (*(*_startFrontEdges[i])->vertexes[0])->belongsToShellVertex &&
			 ((*(*_startFrontEdges[i])->vertexes[1])->belongsToShellEdge)->IsContaining(*(*(*_startFrontEdges[i])->vertexes[0])->belongsToShellVertex)))
			continue;
		

		(*_startFrontEdges[i])->Find2FacetsAround(_startFrontFacets, around_facets[0], around_facets[1]);
		around_nodes[0] = (*around_facets[0])->FindVertexNotIncludedInEdge(**_startFrontEdges[i]);
		around_nodes[1] = (*around_facets[1])->FindVertexNotIncludedInEdge(**_startFrontEdges[i]);
		if ((*around_nodes[0])->belongsToShellEdge && (*around_nodes[1])->belongsToShellVertex ||
			(*around_nodes[1])->belongsToShellEdge && (*around_nodes[0])->belongsToShellVertex)
			continue;

		(*_startFrontEdges[i])->FlipIfNeeded(_startFrontEdges, _startFrontFacets);
	}
}

void Polycrystalline3::TriangulateCrystallitesVolumes()
{
	//#pragma omp parallel for
	for (size_t i = 0, max = crystallites.size(); i < max; i++)
	{
		crystallites[i]->SetStartFront(_startFrontEdges, _startFrontFacets);
		crystallites[i]->TriangulateVolume(_preferredLength);
	}
}

void Polycrystalline3::InputData()
{
	_preferredLength = 0.6;

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


	_startFrontVertexes.insert(_startFrontVertexes.end(),
	{
		(new Vertex3(1.0, 0.0, 0.0))->GetPtrToUniquePtr(),
		(new Vertex3(0.5, 0.0, 1.0))->GetPtrToUniquePtr(),
		(new Vertex3(0.0, 0.0, 0.0))->GetPtrToUniquePtr(),
		(new Vertex3(0.5, 1.0, 0.0))->GetPtrToUniquePtr()
	});
	//for (size_t i = 0, shell_vertexes_num = _shellVertexes.size(); i < shell_vertexes_num; i++)
	//{
	//	(*_startFrontVertexes[i])->belongsToShellVertex = _shellVertexes[i];
	//}
	_startFrontEdges.insert(_startFrontEdges.end(),
	{
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[1]))->GetPtrToUniquePtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[2]))->GetPtrToUniquePtr(),
		(new Edge3(**_startFrontVertexes[2], **_startFrontVertexes[3]))->GetPtrToUniquePtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[3]))->GetPtrToUniquePtr(),
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[3]))->GetPtrToUniquePtr(),
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[2]))->GetPtrToUniquePtr()
	});
	_startFrontFacets.insert(_startFrontFacets.end(),
	{
		(new Facet3(**_startFrontEdges[0], **_startFrontEdges[1], **_startFrontEdges[5]))->GetPtrToUniquePtr(),
		(new Facet3(**_startFrontEdges[1], **_startFrontEdges[2], **_startFrontEdges[3]))->GetPtrToUniquePtr(),
		(new Facet3(**_startFrontEdges[0], **_startFrontEdges[3], **_startFrontEdges[4]))->GetPtrToUniquePtr(),
		(new Facet3(**_startFrontEdges[2], **_startFrontEdges[4], **_startFrontEdges[5]))->GetPtrToUniquePtr()
	});
}

void Polycrystalline3::OutputData(string filename) const
{
	ofstream file(filename);
	
	for (size_t i = 0, verts_num = _startFrontVertexes.size(); i < verts_num; i++)
	{
		file << "v " << (**_startFrontVertexes[i])[0] << ' ' << (**_startFrontVertexes[i])[1] << ' ' << (**_startFrontVertexes[i])[2] << '\n';
		(*_startFrontVertexes[i])->globalNum = i + 1;
	}
	//for (auto &facet : _startFrontFacets)
	//{
	//	if (!*facet)
	//		continue;

	//	vector<size_t> gl_nums;
	//	for (auto &edge : (*facet)->edges)
	//		for (auto &vert : (*edge)->vertexes)
	//			if (std::find(gl_nums.begin(), gl_nums.end(), (*vert)->globalNum) == gl_nums.end())
	//				gl_nums.push_back((*vert)->globalNum);

	//	file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
	//}
	for (auto &f_facet : crystallites[0]->frontFacets)
	{
		if (!*f_facet)
			continue;
		auto facet = (*f_facet)->facet->GetPtrToUniquePtr();

		vector<size_t> gl_nums;
		for (auto &edge : (*facet)->edges)
			for (auto &vert : (*edge)->vertexes)
				if (std::find(gl_nums.begin(), gl_nums.end(), (*vert)->globalNum) == gl_nums.end())
					gl_nums.push_back((*vert)->globalNum);

		file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
	}
	for (auto &facet : crystallites[0]->innerFacets)
	{
		if (!*facet)
			continue;

		vector<size_t> gl_nums;
		for (auto &edge : (*facet)->edges)
			for (auto &vert : (*edge)->vertexes)
				if (std::find(gl_nums.begin(), gl_nums.end(), (*vert)->globalNum) == gl_nums.end())
					gl_nums.push_back((*vert)->globalNum);

		file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
	}

	file.close();
}

Polycrystalline3::Polycrystalline3() {}

Polycrystalline3::~Polycrystalline3()
{
	for (auto &crys : crystallites)
		delete crys;

	for (auto &facet : _shellFacets)
		delete facet;
	for (auto &edge : _shellEdges)
		delete edge;
	for (auto &vert : _shellVertexes)
		delete vert;

	for (auto &ptr : _startFrontFacets)
	{
		delete ptr->release();
		delete ptr;
	}
	for (auto &ptr : _startFrontEdges)
	{
		delete ptr->release();
		delete ptr;
	}
	for (auto &ptr : _startFrontVertexes)
	{
		delete ptr->release();
		delete ptr;
	}
}