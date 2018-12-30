#include "Polycrystal3.h"
#include <algorithm>
#include <iostream>
#include "Timer.h"


ShellEdge3* Polycrystal3::findShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const
{
	for (auto &s_edge : _shellEdges)
		if ((s_edge->vertexes[0] == v0  &&
			 s_edge->vertexes[1] == v1) ||
			(s_edge->vertexes[1] == v0  &&
			 s_edge->vertexes[0] == v1))
		{
			return s_edge;
		}
	
	return nullptr;
}

unique_ptr<Edge3>* Polycrystal3::findStartFrontEdge(const unique_ptr<Vertex3>* v0, const unique_ptr<Vertex3>* v1) const
{
	for (auto &edge : _startFrontEdges)
		if (((*edge)->vertexes[0] == v0  &&
			 (*edge)->vertexes[1] == v1) ||
			((*edge)->vertexes[1] == v0  &&
			 (*edge)->vertexes[0] == v1))
		{
			return edge;
		}

	return nullptr;
}

void Polycrystal3::setLinksWithShell()
{
	double sqr_sufficient_dist = _preferredLength * _preferredLength * 1e-6;
	size_t shell_verts_num = _shellVertexes.size();
	size_t shell_edges_num = _shellEdges.size();
	size_t shell_facets_num = _shellFacets.size();
	#pragma omp parallel firstprivate(sqr_sufficient_dist, shell_verts_num, shell_edges_num, shell_facets_num)
	{
		#pragma omp for
		for (size_t i = 0ull; i < shell_verts_num; i++)
		{
			for (auto &vert : _startFrontVertexes)
			{
				if ((*_shellVertexes[i] - **vert).sqrMagnitude() < sqr_sufficient_dist)
				{
					(*vert)->belongsToShellVertex = _shellVertexes[i];
					break;
				}
			}
		}
		#pragma omp for
		for (size_t i = 0ull; i < shell_edges_num; i++)
		{
			Point3 proj_buf;
			for (auto &vert : _startFrontVertexes)
			{
				if ((*vert)->belongsToShellVertex)
					continue;

				if (spatialalgs::project(
						proj_buf,
						(*vert)->getPosition(),
						_shellEdges[i]->vertexes[0]->getPosition(), _shellEdges[i]->vertexes[1]->getPosition()) &&
						(proj_buf - (*vert)->getPosition()).sqrMagnitude() < sqr_sufficient_dist)
				{
					(*vert)->belongsToShellEdge = _shellEdges[i];
				}
			}
		}
		#pragma omp for
		for (size_t i = 0ull; i < shell_facets_num; i++)
		{
			Point3 proj_buf;
			for (auto &vert : _startFrontVertexes)
			{
				if ((*vert)->belongsToShellVertex ||
					(*vert)->belongsToShellEdge)
					continue;

				proj_buf = _shellFacets[i]->edges[0]->vertexes[0]->getPosition()
							+ (**vert - *_shellFacets[i]->edges[0]->vertexes[0]).project(
									*_shellFacets[i]->edges[0]->vertexes[1] - *_shellFacets[i]->edges[0]->vertexes[0],
									*_shellFacets[i]->edges[1]->vertexes[1] - *_shellFacets[i]->edges[1]->vertexes[0]);

				if ((proj_buf - (*vert)->getPosition()).sqrMagnitude() < sqr_sufficient_dist)
					(*vert)->belongsToShellFacet = _shellFacets[i];
			}
		}
	}
}

template <class T>
void Polycrystal3::removePtrsToNullptr(vector<unique_ptr<T>*> &vec)
{
	size_t real_objs_num = std::count_if(vec.begin(), vec.end(), [](unique_ptr<T>*& ptr) { return *ptr ? true : false; });
	vector<unique_ptr<T>*> buf_vec(real_objs_num);

	size_t firts_thr_nums = real_objs_num / 2ull;
	size_t second_thr_nums = real_objs_num - firts_thr_nums;
	//#pragma omp parallel num_threads(2)
	//{
	//#pragma omp single
	//{
	size_t index1 = 0ull;
	for (size_t i = 0ull; index1 < firts_thr_nums; i++)
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
	size_t index2 = 0ull;
	for (size_t i = vec.size() - 1ull; index2 < second_thr_nums; i--)
	{
		if (*vec[i])
		{
			buf_vec[real_objs_num - 1ull - index2] = vec[i];
			index2++;
		}
	}
	//}
	//}

	vec = std::move(buf_vec);
}

void Polycrystal3::removePtrsToNullptrFromVectors()
{
	//#pragma omp parallel num_threads(3)
	//{
	//#pragma omp single
	removePtrsToNullptr(_startFrontVertexes);
	//#pragma omp single
	removePtrsToNullptr(_startFrontEdges);
	//#pragma omp single
	removePtrsToNullptr(_startFrontFacets);
	//}
}

void Polycrystal3::triangulateShell()
{
	size_t edges_num = _startFrontEdges.size();
	for (size_t i = 0ull; i < edges_num; i++)
	{
		if (*_startFrontEdges[i] &&
			(*_startFrontEdges[i])->sqrMagnitude() > 2.25 * _preferredLength * _preferredLength)
		{
			(*_startFrontEdges[i])->make2Instead(_startFrontFacets, _startFrontEdges, _startFrontVertexes);
			edges_num += 2ull;
		}
	}

	removePtrsToNullptrFromVectors();
}

void Polycrystal3::startFrontDelaunayPostprocessing()
{
	unique_ptr<Vertex3>* around_nodes[2];
	unique_ptr<Facet3>* around_facets[2];
	size_t edges_num = _startFrontEdges.size();
	for (size_t i = 0ull; i < edges_num; i++)
	{
		if (!*_startFrontEdges[i] ||
			((*(*_startFrontEdges[i])->vertexes[0])->belongsToShellEdge == (*(*_startFrontEdges[i])->vertexes[1])->belongsToShellEdge && 
			 (*(*_startFrontEdges[i])->vertexes[1])->belongsToShellEdge) ||
			((*(*_startFrontEdges[i])->vertexes[0])->belongsToShellEdge && 
			 (*(*_startFrontEdges[i])->vertexes[1])->belongsToShellVertex && 
			 (*(*_startFrontEdges[i])->vertexes[0])->belongsToShellEdge->contains(*(*(*_startFrontEdges[i])->vertexes[1])->belongsToShellVertex)) ||
			((*(*_startFrontEdges[i])->vertexes[1])->belongsToShellEdge &&
			 (*(*_startFrontEdges[i])->vertexes[0])->belongsToShellVertex &&
			 ((*(*_startFrontEdges[i])->vertexes[1])->belongsToShellEdge)->contains(*(*(*_startFrontEdges[i])->vertexes[0])->belongsToShellVertex)))
			continue;
		

		(*_startFrontEdges[i])->find2AdjFacets(_startFrontFacets, around_facets[0], around_facets[1]);
		around_nodes[0] = (*around_facets[0])->findVertexNotIncludedInEdge(**_startFrontEdges[i]);
		around_nodes[1] = (*around_facets[1])->findVertexNotIncludedInEdge(**_startFrontEdges[i]);
		if ((*around_nodes[0])->belongsToShellEdge && (*around_nodes[1])->belongsToShellVertex ||
			(*around_nodes[1])->belongsToShellEdge && (*around_nodes[0])->belongsToShellVertex)
			continue;

		(*_startFrontEdges[i])->flipIfNeeded(_startFrontEdges, _startFrontFacets);
	}
}

PolycrMesh* Polycrystal3::structurizeMesh()
{
	PolycrMesh* pol_triang = new PolycrMesh;

	pol_triang->crysesNum = crystallites.size();
	pol_triang->crysesTetrsNum = new size_t[pol_triang->crysesNum];

	size_t nodes_num = _startFrontVertexes.size();
	for (auto &crys : crystallites)
		nodes_num += std::count_if(
			crys->innerVerts.begin(),
			crys->innerVerts.end(),
			[](auto vert)-> bool { return (bool)*vert; });

	pol_triang->nodesNum = nodes_num;
	pol_triang->nodesPositions = new double[3ull * pol_triang->nodesNum];
	unique_ptr<Vertex3>** nodes_ptrs = new unique_ptr<Vertex3>*[pol_triang->nodesNum];
	auto FindNodeIndex = [&pol_triang, nodes_ptrs](unique_ptr<Vertex3>* node_ptr)-> size_t
	{
		for (size_t i = 0ull; i < pol_triang->nodesNum; i++)
		{
			if (nodes_ptrs[i] == node_ptr)
				return i;
		}
		return pol_triang->nodesNum;
	};

	size_t node_ind = 0ull;
	for (auto &s_f_vert : _startFrontVertexes)
	{
		if (!*s_f_vert)
			continue;

		nodes_ptrs[node_ind] = s_f_vert;
		pol_triang->nodesPositions[3ull * node_ind]        = (**s_f_vert)[0];
		pol_triang->nodesPositions[3ull * node_ind + 1ull] = (**s_f_vert)[1];
		pol_triang->nodesPositions[3ull * node_ind + 2ull] = (**s_f_vert)[2];
		node_ind++;
	}
	for (auto &crys : crystallites)
	{
		for (auto &vert : crys->innerVerts)
		{
			if (!*vert)
				continue;

			nodes_ptrs[node_ind] = vert;
			pol_triang->nodesPositions[3ull * node_ind]        = (**vert)[0];
			pol_triang->nodesPositions[3ull * node_ind + 1ull] = (**vert)[1];
			pol_triang->nodesPositions[3ull * node_ind + 2ull] = (**vert)[2];
			node_ind++;
		}
	}

	size_t tetrs_num = 0ull;
	for (size_t i = 0ull; i < pol_triang->crysesNum; i++)
	{
		pol_triang->crysesTetrsNum[i] = std::count_if(
			crystallites[i]->innerSimps.begin(),
			crystallites[i]->innerSimps.end(),
			[](auto simp)-> bool { return (bool)*simp; });

		tetrs_num += pol_triang->crysesTetrsNum[i];
	}

	pol_triang->tetrsNum = tetrs_num;
	pol_triang->tetrs = new size_t[4ull * pol_triang->tetrsNum];

	size_t tetr_ind = 0ull;
	for (size_t i = 0ull; i < pol_triang->crysesNum; i++)
	{
		for (auto &tetr : crystallites[i]->innerSimps)
		{
			if (!*tetr)
				continue;

			pol_triang->tetrs[4ull * tetr_ind]        = FindNodeIndex((*tetr)->vertexes[0]);
			pol_triang->tetrs[4ull * tetr_ind + 1ull] = FindNodeIndex((*tetr)->vertexes[1]);
			pol_triang->tetrs[4ull * tetr_ind + 2ull] = FindNodeIndex((*tetr)->vertexes[1]);
			pol_triang->tetrs[4ull * tetr_ind + 3ull] = FindNodeIndex((*tetr)->vertexes[1]);
			tetr_ind++;
		}
	}

	return pol_triang;
}

void Polycrystal3::generateMeshNoStruct(const double preferredLength)
{
	_preferredLength = preferredLength;

	triangulateShell();
	setLinksWithShell();
	startFrontDelaunayPostprocessing();

	Timer tmr;
	tmr.Start();
	double min_q = 0.0, av_q = 0.0;
	#pragma omp parallel for
	for (size_t i = 0ull, max = crystallites.size(); i < max; i++)
	{
		crystallites[i]->setStartFront(_startFrontEdges, _startFrontFacets);
		crystallites[i]->generateMesh(_preferredLength, this);

		double buf_min_q, buf_av_q;
		crystallites[i]->analyzeMeshQuality(buf_min_q, buf_av_q);
		min_q += buf_min_q;
		av_q += buf_av_q;
	}
	min_q /= crystallites.size();
	av_q /= crystallites.size();
	tmr.Stop();
	std::ofstream out("log.txt");
	out << ">>  Minimum quality    " << min_q
		<< "\n>>  Average quality    " << av_q
		<< "\n>>  Run time           " << tmr.GetLastDuration(microseconds) * 1e-6 << " s";
	out.close();
}

void Polycrystal3::generateMeshNoStruct(string filename, const double preferredLength)
{
	inputData(filename);
	generateMeshNoStruct(preferredLength);
}

void Polycrystal3::generateMeshNoStruct(const CrysesShell& crysesShell, const double preferredLength)
{
	inputData(crysesShell);
	generateMeshNoStruct(preferredLength);
}

PolycrMesh* Polycrystal3::generateMesh(const double preferredLength)
{
	generateMeshNoStruct(preferredLength);
	return structurizeMesh();
}

PolycrMesh* Polycrystal3::generateMesh(string filename, const double preferredLength)
{
	generateMeshNoStruct(filename, preferredLength);
	return structurizeMesh();
}

PolycrMesh* Polycrystal3::generateMesh(const CrysesShell& crysesShell, const double preferredLength)
{
	generateMeshNoStruct(crysesShell, preferredLength);
	return structurizeMesh();
}

PolycrMesh* Polycrystal3::getLastMesh()
{
	return lastTriangulation;
}

void Polycrystal3::inputData()
{
	_shellVertexes.insert(_shellVertexes.end(),
	{
		new ShellVertex3(1.0, 0.0, 0.0),
		new ShellVertex3(1.0, 1.0, 0.0),
		new ShellVertex3(0.0, 0.0, 0.0),
		new ShellVertex3(1.0, 0.0, 1.0),
		new ShellVertex3(0.0, 0.0, 1.0),
		new ShellVertex3(1.0, 1.0, 1.0),
		new ShellVertex3(0.0, 1.0, 0.0),
		new ShellVertex3(0.0, 1.0, 1.0)
	});
	_shellEdges.insert(_shellEdges.end(),
	{
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[1]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[2]),
		new ShellEdge3(*_shellVertexes[2], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[2]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[5]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[6]),
		new ShellEdge3(*_shellVertexes[2], *_shellVertexes[4]),
		new ShellEdge3(*_shellVertexes[3], *_shellVertexes[4]),
		new ShellEdge3(*_shellVertexes[4], *_shellVertexes[5]),
		new ShellEdge3(*_shellVertexes[3], *_shellVertexes[5]),
		new ShellEdge3(*_shellVertexes[2], *_shellVertexes[6]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[7]),
		new ShellEdge3(*_shellVertexes[6], *_shellVertexes[7]),
		new ShellEdge3(*_shellVertexes[5], *_shellVertexes[7]),
		new ShellEdge3(*_shellVertexes[2], *_shellVertexes[7]),
		new ShellEdge3(*_shellVertexes[4], *_shellVertexes[7])
	});
	_shellFacets.insert(_shellFacets.end(),
	{
		new ShellFacet3(*_shellEdges[0], *_shellEdges[1], *_shellEdges[5]),
		new ShellFacet3(*_shellEdges[0], *_shellEdges[3], *_shellEdges[4]),
		new ShellFacet3(*_shellEdges[2], *_shellEdges[4], *_shellEdges[5]),
		new ShellFacet3(*_shellEdges[3], *_shellEdges[6], *_shellEdges[11]),
		new ShellFacet3(*_shellEdges[2], *_shellEdges[8], *_shellEdges[9]),
		new ShellFacet3(*_shellEdges[6], *_shellEdges[13], *_shellEdges[15]),
		new ShellFacet3(*_shellEdges[10], *_shellEdges[15], *_shellEdges[17]),
		new ShellFacet3(*_shellEdges[9], *_shellEdges[10], *_shellEdges[11]),
		new ShellFacet3(*_shellEdges[7], *_shellEdges[13], *_shellEdges[14]),
		new ShellFacet3(*_shellEdges[1], *_shellEdges[7], *_shellEdges[12]),
		new ShellFacet3(*_shellEdges[12], *_shellEdges[14], *_shellEdges[16]),
		new ShellFacet3(*_shellEdges[8], *_shellEdges[16], *_shellEdges[17])
	});

	crystallites.push_back(new Crystallite3);
	crystallites.front()->shellEdges = _shellEdges;
	crystallites.front()->shellFacets = _shellFacets;


	_startFrontVertexes.insert(_startFrontVertexes.end(),
	{
		(new Vertex3(1.0, 0.0, 0.0))->getPtrToUPtr(),
		(new Vertex3(1.0, 1.0, 0.0))->getPtrToUPtr(),
		(new Vertex3(0.0, 0.0, 0.0))->getPtrToUPtr(),
		(new Vertex3(1.0, 0.0, 1.0))->getPtrToUPtr(),
		(new Vertex3(0.0, 0.0, 1.0))->getPtrToUPtr(),
		(new Vertex3(1.0, 1.0, 1.0))->getPtrToUPtr(),
		(new Vertex3(0.0, 1.0, 0.0))->getPtrToUPtr(),
		(new Vertex3(0.0, 1.0, 1.0))->getPtrToUPtr()
	});
	_startFrontEdges.insert(_startFrontEdges.end(),
	{
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[1]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[2]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[2], **_startFrontVertexes[3]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[3]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[3]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[2]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[5]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[6]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[2], **_startFrontVertexes[4]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[3], **_startFrontVertexes[4]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[4], **_startFrontVertexes[5]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[3], **_startFrontVertexes[5]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[2], **_startFrontVertexes[6]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[7]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[6], **_startFrontVertexes[7]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[5], **_startFrontVertexes[7]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[2], **_startFrontVertexes[7]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[4], **_startFrontVertexes[7]))->getPtrToUPtr()
	});
	_startFrontFacets.insert(_startFrontFacets.end(),
	{
		(new Facet3(**_startFrontEdges[0], **_startFrontEdges[1], **_startFrontEdges[5]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[0], **_startFrontEdges[3], **_startFrontEdges[4]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[2], **_startFrontEdges[4], **_startFrontEdges[5]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[3], **_startFrontEdges[6], **_startFrontEdges[11]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[2], **_startFrontEdges[8], **_startFrontEdges[9]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[6], **_startFrontEdges[13], **_startFrontEdges[15]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[10], **_startFrontEdges[15], **_startFrontEdges[17]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[9], **_startFrontEdges[10], **_startFrontEdges[11]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[7], **_startFrontEdges[13], **_startFrontEdges[14]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[1], **_startFrontEdges[7], **_startFrontEdges[12]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[12], **_startFrontEdges[14], **_startFrontEdges[16]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[8], **_startFrontEdges[16], **_startFrontEdges[17]))->getPtrToUPtr()
	});

	/*_shellVertexes.insert(_shellVertexes.end(), 
	{
		new ShellVertex3(1.0, 0.0, 0.0),
		new ShellVertex3(0.5, 1.0, 0.0),
		new ShellVertex3(0.0, 0.0, 0.0),
		new ShellVertex3(1.0, 0.0, 1.0),
		new ShellVertex3(0.0, 0.0, 1.0),
		new ShellVertex3(0.5, 1.0, 1.0)
	});
	_shellEdges.insert(_shellEdges.end(),
	{
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[1]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[2]),
		new ShellEdge3(*_shellVertexes[2], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[3]),
		new ShellEdge3(*_shellVertexes[0], *_shellVertexes[2]),
		new ShellEdge3(*_shellVertexes[1], *_shellVertexes[5]),
		new ShellEdge3(*_shellVertexes[2], *_shellVertexes[5]),
		new ShellEdge3(*_shellVertexes[2], *_shellVertexes[4]),
		new ShellEdge3(*_shellVertexes[3], *_shellVertexes[4]),
		new ShellEdge3(*_shellVertexes[4], *_shellVertexes[5]),
		new ShellEdge3(*_shellVertexes[3], *_shellVertexes[5])
	});
	_shellFacets.insert(_shellFacets.end(),
	{
		new ShellFacet3(*_shellEdges[0], *_shellEdges[1], *_shellEdges[5]),
		new ShellFacet3(*_shellEdges[0], *_shellEdges[3], *_shellEdges[4]),
		new ShellFacet3(*_shellEdges[2], *_shellEdges[4], *_shellEdges[5]),
		new ShellFacet3(*_shellEdges[3], *_shellEdges[6], *_shellEdges[11]),
		new ShellFacet3(*_shellEdges[2], *_shellEdges[8], *_shellEdges[9]),
		new ShellFacet3(*_shellEdges[1], *_shellEdges[6], *_shellEdges[7]),
		new ShellFacet3(*_shellEdges[7], *_shellEdges[8], *_shellEdges[10]),
		new ShellFacet3(*_shellEdges[9], *_shellEdges[10], *_shellEdges[11])
	});

	crystallites.push_back(new Crystallite3);
	crystallites.front()->shellEdges = _shellEdges;
	crystallites.front()->shellFacets = _shellFacets;


	_startFrontVertexes.insert(_startFrontVertexes.end(),
	{
		(new Vertex3(1.0, 0.0, 0.0))->getPtrToUPtr(),
		(new Vertex3(0.5, 1.0, 0.0))->getPtrToUPtr(),
		(new Vertex3(0.0, 0.0, 0.0))->getPtrToUPtr(),
		(new Vertex3(1.0, 0.0, 1.0))->getPtrToUPtr(),
		(new Vertex3(0.0, 0.0, 1.0))->getPtrToUPtr(),
		(new Vertex3(0.5, 1.0, 1.0))->getPtrToUPtr()
	});
	_startFrontEdges.insert(_startFrontEdges.end(),
	{
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[1]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[2]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[2], **_startFrontVertexes[3]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[3]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[3]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[0], **_startFrontVertexes[2]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[1], **_startFrontVertexes[5]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[2], **_startFrontVertexes[5]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[2], **_startFrontVertexes[4]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[3], **_startFrontVertexes[4]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[4], **_startFrontVertexes[5]))->getPtrToUPtr(),
		(new Edge3(**_startFrontVertexes[3], **_startFrontVertexes[5]))->getPtrToUPtr()
	});
	_startFrontFacets.insert(_startFrontFacets.end(),
	{
		(new Facet3(**_startFrontEdges[0], **_startFrontEdges[1], **_startFrontEdges[5]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[0], **_startFrontEdges[3], **_startFrontEdges[4]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[2], **_startFrontEdges[4], **_startFrontEdges[5]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[3], **_startFrontEdges[6], **_startFrontEdges[11]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[2], **_startFrontEdges[8], **_startFrontEdges[9]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[1], **_startFrontEdges[6], **_startFrontEdges[7]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[7], **_startFrontEdges[8], **_startFrontEdges[10]))->getPtrToUPtr(),
		(new Facet3(**_startFrontEdges[9], **_startFrontEdges[10], **_startFrontEdges[11]))->getPtrToUPtr()
	});*/
}

void Polycrystal3::inputData(string filename)
{
	ifstream input(filename);

	size_t nodes_num;
	input >> nodes_num;

	for (size_t i = 0ull; i < nodes_num; i++)
	{
		double coors[3];
		input >> coors[0];
		input >> coors[1];
		input >> coors[2];

		_shellVertexes.push_back(new ShellVertex3(coors[0], coors[1], coors[2]));
		_startFrontVertexes.push_back((new Vertex3(coors[0], coors[1], coors[2]))->getPtrToUPtr());
	}

	size_t facets_num;
	input >> facets_num;

	for (size_t i = 0ull; i < facets_num; i++)
	{
		size_t facet_nodes_inds[3];
		input >> facet_nodes_inds[0];
		input >> facet_nodes_inds[1];
		input >> facet_nodes_inds[2];

		if (!findShellEdge(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]))
		{
			_shellEdges.push_back(new ShellEdge3(*_shellVertexes[facet_nodes_inds[0]], *_shellVertexes[facet_nodes_inds[1]]));
			_startFrontEdges.push_back((new Edge3(**_startFrontVertexes[facet_nodes_inds[0]], **_startFrontVertexes[facet_nodes_inds[1]]))->getPtrToUPtr());
		}

		if (!findShellEdge(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]))
		{
			_shellEdges.push_back(new ShellEdge3(*_shellVertexes[facet_nodes_inds[1]], *_shellVertexes[facet_nodes_inds[2]]));
			_startFrontEdges.push_back((new Edge3(**_startFrontVertexes[facet_nodes_inds[1]], **_startFrontVertexes[facet_nodes_inds[2]]))->getPtrToUPtr());
		}

		if (!findShellEdge(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]]))
		{
			_shellEdges.push_back(new ShellEdge3(*_shellVertexes[facet_nodes_inds[2]], *_shellVertexes[facet_nodes_inds[0]]));
			_startFrontEdges.push_back((new Edge3(**_startFrontVertexes[facet_nodes_inds[2]], **_startFrontVertexes[facet_nodes_inds[0]]))->getPtrToUPtr());
		}

		_shellFacets.push_back(new ShellFacet3(
			*findShellEdge(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]),
			*findShellEdge(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]),
			*findShellEdge(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]])));
		_startFrontFacets.push_back((new Facet3(
			**findStartFrontEdge(_startFrontVertexes[facet_nodes_inds[0]], _startFrontVertexes[facet_nodes_inds[1]]),
			**findStartFrontEdge(_startFrontVertexes[facet_nodes_inds[1]], _startFrontVertexes[facet_nodes_inds[2]]),
			**findStartFrontEdge(_startFrontVertexes[facet_nodes_inds[2]], _startFrontVertexes[facet_nodes_inds[0]])))->getPtrToUPtr());
	}

	size_t cryses_num;
	input >> cryses_num;
	crystallites.insert(crystallites.end(), cryses_num, new Crystallite3);

	size_t* cryses_facets_nums = new size_t[cryses_num];
	for (size_t i = 0ull; i < cryses_num; i++)
		input >> cryses_facets_nums[i];

	for (size_t i = 0ull; i < cryses_num; i++)
	{
		for (size_t j = 0ull; j < cryses_facets_nums[i]; j++)
		{
			size_t facet_ind;
			input >> facet_ind;

			if (std::find(crystallites[i]->shellEdges.begin(), crystallites[i]->shellEdges.end(), _shellFacets[facet_ind]->edges[0]) == crystallites[i]->shellEdges.end())
				crystallites[i]->shellEdges.push_back(_shellFacets[facet_ind]->edges[0]);

			if (std::find(crystallites[i]->shellEdges.begin(), crystallites[i]->shellEdges.end(), _shellFacets[facet_ind]->edges[1]) == crystallites[i]->shellEdges.end())
				crystallites[i]->shellEdges.push_back(_shellFacets[facet_ind]->edges[1]);

			if (std::find(crystallites[i]->shellEdges.begin(), crystallites[i]->shellEdges.end(), _shellFacets[facet_ind]->edges[2]) == crystallites[i]->shellEdges.end())
				crystallites[i]->shellEdges.push_back(_shellFacets[facet_ind]->edges[2]);

			crystallites[i]->shellFacets.push_back(_shellFacets[facet_ind]);
		}
	}

	input.close();
	delete[] cryses_facets_nums;
}

void Polycrystal3::inputData(const CrysesShell& crysesShell)
{
	for (size_t i = 0ull; i < crysesShell.nodesNum; i++)
	{
		double coors[3];
		coors[0] = crysesShell.nodesPositions[3ull * i];
		coors[1] = crysesShell.nodesPositions[3ull * i + 1ull];
		coors[2] = crysesShell.nodesPositions[3ull * i + 2ull];

		_shellVertexes.push_back(new ShellVertex3(coors[0], coors[1], coors[2]));
		_startFrontVertexes.push_back((new Vertex3(coors[0], coors[1], coors[2]))->getPtrToUPtr());
	}

	for (size_t i = 0ull; i < crysesShell.facetsNum; i++)
	{
		size_t facet_nodes_inds[3];
		facet_nodes_inds[0] = crysesShell.facets[3ull * i];
		facet_nodes_inds[1] = crysesShell.facets[3ull * i + 1ull];
		facet_nodes_inds[2] = crysesShell.facets[3ull * i + 2ull];

		if (!findShellEdge(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]))
		{
			_shellEdges.push_back(new ShellEdge3(*_shellVertexes[facet_nodes_inds[0]], *_shellVertexes[facet_nodes_inds[1]]));
			_startFrontEdges.push_back((new Edge3(**_startFrontVertexes[facet_nodes_inds[0]], **_startFrontVertexes[facet_nodes_inds[1]]))->getPtrToUPtr());
		}

		if (!findShellEdge(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]))
		{
			_shellEdges.push_back(new ShellEdge3(*_shellVertexes[facet_nodes_inds[1]], *_shellVertexes[facet_nodes_inds[2]]));
			_startFrontEdges.push_back((new Edge3(**_startFrontVertexes[facet_nodes_inds[1]], **_startFrontVertexes[facet_nodes_inds[2]]))->getPtrToUPtr());
		}

		if (!findShellEdge(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]]))
		{
			_shellEdges.push_back(new ShellEdge3(*_shellVertexes[facet_nodes_inds[2]], *_shellVertexes[facet_nodes_inds[0]]));
			_startFrontEdges.push_back((new Edge3(**_startFrontVertexes[facet_nodes_inds[2]], **_startFrontVertexes[facet_nodes_inds[0]]))->getPtrToUPtr());
		}

		_shellFacets.push_back(new ShellFacet3(
			*findShellEdge(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]),
			*findShellEdge(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]),
			*findShellEdge(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]])));
		_startFrontFacets.push_back((new Facet3(
			**findStartFrontEdge(_startFrontVertexes[facet_nodes_inds[0]], _startFrontVertexes[facet_nodes_inds[1]]),
			**findStartFrontEdge(_startFrontVertexes[facet_nodes_inds[1]], _startFrontVertexes[facet_nodes_inds[2]]),
			**findStartFrontEdge(_startFrontVertexes[facet_nodes_inds[2]], _startFrontVertexes[facet_nodes_inds[0]])))->getPtrToUPtr());
	}

	crystallites.insert(crystallites.end(), crysesShell.crysesNum, new Crystallite3);
	
	size_t cryses_inner_ind = 0ull;
	for (size_t i = 0ull; i < crysesShell.crysesNum; i++)
	{
		for (size_t j = 0ull; j < crysesShell.crysesFacetsNums[i]; j++)
		{
			size_t facet_ind = crysesShell.cryses[cryses_inner_ind++];

			if (std::find(crystallites[i]->shellEdges.begin(), crystallites[i]->shellEdges.end(), _shellFacets[facet_ind]->edges[0]) == crystallites[i]->shellEdges.end())
				crystallites[i]->shellEdges.push_back(_shellFacets[facet_ind]->edges[0]);

			if (std::find(crystallites[i]->shellEdges.begin(), crystallites[i]->shellEdges.end(), _shellFacets[facet_ind]->edges[1]) == crystallites[i]->shellEdges.end())
				crystallites[i]->shellEdges.push_back(_shellFacets[facet_ind]->edges[1]);

			if (std::find(crystallites[i]->shellEdges.begin(), crystallites[i]->shellEdges.end(), _shellFacets[facet_ind]->edges[2]) == crystallites[i]->shellEdges.end())
				crystallites[i]->shellEdges.push_back(_shellFacets[facet_ind]->edges[2]);

			crystallites[i]->shellFacets.push_back(_shellFacets[facet_ind]);
		}
	}
}

void Polycrystal3::outputData(string filename) const
{
	ofstream file(filename);
	size_t i = 0ull;
	for (size_t verts_num = _startFrontVertexes.size(); i < verts_num; i++)
	{
		file << "v " << (**_startFrontVertexes[i])[0] << ' ' << (**_startFrontVertexes[i])[1] << ' ' << (**_startFrontVertexes[i])[2] << '\n';
		(*_startFrontVertexes[i])->globalNum = i + 1ull;
	}
	for (auto &crys : crystallites)
	{
		if (!crys)
			continue;

		for (size_t j = 0ull, verts_num = crys->innerVerts.size(); j < verts_num; i++, j++)
		{
			file << "v " << (**crys->innerVerts[j])[0] << ' ' << (**crys->innerVerts[j])[1] << ' ' << (**crys->innerVerts[j])[2] << '\n';
			(*crys->innerVerts[j])->globalNum = i + 1ull;
		}
	}
#define DEV_DEBUG
#ifdef DEV_DEBUG
	for (auto &f_facet : crystallites[0]->frontFacets)
	{
		if (!*f_facet)
			continue;
		auto facet = (*f_facet)->facet->getPtrToUPtr();
	
		vector<size_t> gl_nums;
		for (auto &edge : (*facet)->edges)
			for (auto &vert : (*edge)->vertexes)
				if (std::find(gl_nums.begin(), gl_nums.end(), (*vert)->globalNum) == gl_nums.end())
					gl_nums.push_back((*vert)->globalNum);
	
		file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
	}
#else
	for (auto &facet : _startFrontFacets)
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
#endif

	file.close();
}

Polycrystal3::Polycrystal3() {}

Polycrystal3::Polycrystal3(string filename)
{
	inputData(filename);
}

Polycrystal3::Polycrystal3(const CrysesShell& crysesShell)
{
	inputData(crysesShell);
}

Polycrystal3::~Polycrystal3()
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