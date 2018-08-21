#include "Polycrystalline2.h"
#include <algorithm>
#include <iostream>

#define DET(a, b, c, d) \
		(a * d - b * c)

#define EPS 1e-10
#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - EPS < p && p < std::max(p0_coor, p1_coor) + EPS)

#define INSIDE_RECTANGLE(corner0, corner1, point) \
		(BETWEEN(corner0[0], corner1[0], point[0]) && \
		 BETWEEN(corner0[1], corner1[1], point[1]))

#define NOT_INTERSECT_2ND_CHECK(edge0node0, edge0node1, edge1node0, edge1node1) \
		(edge0node0[0] < edge1node0[0]  && \
		 edge0node1[0] < edge1node0[0] && \
		 edge0node0[0] < edge1node1[0] && \
		 edge0node1[0] < edge1node1[0]) || \
		\
		(edge0node0[1] < edge1node0[1] && \
		 edge0node1[1] < edge1node0[1] && \
		 edge0node0[1] < edge1node1[1] && \
		 edge0node1[1] < edge1node1[1]) || \
		\
		(edge0node0[0] > edge1node0[0]  && \
		 edge0node1[0] > edge1node0[0] && \
		 edge0node0[0] > edge1node1[0] && \
		 edge0node1[0] > edge1node1[0]) || \
		\
		(edge0node0[1] > edge1node0[1] && \
		 edge0node1[1] > edge1node0[1] && \
		 edge0node0[1] > edge1node1[1] && \
		 edge0node1[1] > edge1node1[1])

void Polycrystalline2::Debug()
{
}

void Polycrystalline2::GenerateNodesEvenly(double* const polycrysSizeAxis, size_t* const minNodesNumAxis)
{
	double startGenCoor[2];
	startGenCoor[0] = _minShellNodesCoor[0] + (polycrysSizeAxis[0] - minNodesNumAxis[0] * _l_av) * 0.5;
	startGenCoor[1] = _minShellNodesCoor[1] + (polycrysSizeAxis[1] - (minNodesNumAxis[1] - 1) * 0.5 * sqrt(3) * _l_av) * 0.5;

	size_t freeNodesNum = minNodesNumAxis[0] * minNodesNumAxis[1] + (minNodesNumAxis[1] + 1) / 2;

	double node_coors[2];
	node_coors[0] = startGenCoor[0];
	node_coors[1] = startGenCoor[1];
	double delta_node_coors_1 = 0.5 * sqrt(3) * _l_av;
	for (size_t i = 0; i < freeNodesNum;)
	{
		_freeNodes.push_back(
			(new Node2(
				node_coors[0],
				node_coors[1]))
			->GetPtrToUniquePtr());
		i++;
		for (size_t j = 0; j < minNodesNumAxis[0]; j++)
		{
			node_coors[0] += _l_av;
			_freeNodes.push_back(
				(new Node2(
					node_coors[0],
					node_coors[1]))
				->GetPtrToUniquePtr());
			i++;
		}

		if (!(i < freeNodesNum))
		{
			break;
		}

		node_coors[0] -= minNodesNumAxis[0] * _l_av - _l_av * 0.5;
		node_coors[1] += delta_node_coors_1;
		_freeNodes.push_back(
			(new Node2(
				node_coors[0],
				node_coors[1]))
			->GetPtrToUniquePtr());
		i++;
		for (size_t j = 0; j < minNodesNumAxis[0] - 1; j++)
		{
			node_coors[0] += _l_av;
			_freeNodes.push_back(
				(new Node2(
					node_coors[0],
					node_coors[1]))
				->GetPtrToUniquePtr());
			i++;
		}
		node_coors[0] -= minNodesNumAxis[0] * _l_av - _l_av * 0.5;
		node_coors[1] += delta_node_coors_1;
	}
}

unique_ptr<Edge2>* Polycrystalline2::FindEdge(Node2& node0, Node2& node1)
{
	if (std::find(node0.neighbors.begin(), node0.neighbors.end(), node1.GetPtrToUniquePtr()) == node0.neighbors.end())
	{
		return nullptr;
	}

	for (list<unique_ptr<Edge2>*>::const_iterator edges_iter = node0.inclInEdges.begin(), end = node0.inclInEdges.end();
		edges_iter != end; 
		edges_iter++)
	{
		if ((**edges_iter)->nodes[0]->get() == &node0 && 
			(**edges_iter)->nodes[1]->get() == &node1 ||
			(**edges_iter)->nodes[0]->get() == &node1 && 
			(**edges_iter)->nodes[1]->get() == &node0)
		{
			return *edges_iter;
		}
	}

	throw std::exception("Error in function Polycrystalline2::FindFreeEdge");
	return nullptr;
}

void Polycrystalline2::GenerateSimplexesFromNodes(size_t minNodesNumAxis_0)
{
	unique_ptr<Simplex2>* simp_buf = nullptr;
	
	for (size_t i = 0, max = _freeNodes.size() - minNodesNumAxis_0 - 1; i < max;)
	{
		for (size_t j = 0; j < minNodesNumAxis_0; j++)
		{
			_freeSimplexes.push_back(simp_buf = (new Simplex2())->GetPtrToUniquePtr());

			(*simp_buf)->edges[0] = 
				(new Edge2(
					**_freeNodes[i], 
					**_freeNodes[i + 1]))
				->GetPtrToUniquePtr();

			(*simp_buf)->edges[1] = 
				(new Edge2(
					**_freeNodes[i + 1], 
					**_freeNodes[i + 1 + minNodesNumAxis_0]))
				->GetPtrToUniquePtr();

			(*simp_buf)->edges[2] = 
				(new Edge2(
					**_freeNodes[i + 1 + minNodesNumAxis_0], 
					**_freeNodes[i]))
				->GetPtrToUniquePtr();

			_freeEdges.insert(
				_freeEdges.end(),
				{ (*simp_buf)->edges[0],
				  (*simp_buf)->edges[1],
				  (*simp_buf)->edges[2] });

			(*(*simp_buf)->edges[0])->inclInSimplexes.push_back(simp_buf);
			(*(*simp_buf)->edges[1])->inclInSimplexes.push_back(simp_buf);
			(*(*simp_buf)->edges[2])->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;

		if (!(i < max))
		{
			break;
		}

		for (size_t j = 0; j < minNodesNumAxis_0 - 1; j++)
		{
			_freeSimplexes.push_back(simp_buf = (new Simplex2())->GetPtrToUniquePtr());

			(*simp_buf)->edges[0] = 
				(new Edge2(
					**_freeNodes[i], 
					**_freeNodes[i + 1]))
				->GetPtrToUniquePtr();

			(*simp_buf)->edges[1] = 
				(new Edge2(
					**_freeNodes[i + 1], 
					**_freeNodes[i + 1 + minNodesNumAxis_0]))
				->GetPtrToUniquePtr();

			(*simp_buf)->edges[2] = 
				(new Edge2(
					**_freeNodes[i + 1 + minNodesNumAxis_0], 
					**_freeNodes[i]))
				->GetPtrToUniquePtr();

			_freeEdges.insert(
				_freeEdges.end(), 
				{ (*simp_buf)->edges[0],
				  (*simp_buf)->edges[1],
				  (*simp_buf)->edges[2] });

			(*(*simp_buf)->edges[0])->inclInSimplexes.push_back(simp_buf);
			(*(*simp_buf)->edges[1])->inclInSimplexes.push_back(simp_buf);
			(*(*simp_buf)->edges[2])->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;
	}

	unique_ptr<Edge2>* edge_buf = nullptr;
	
	for (size_t i = 1, max = _freeNodes.size() - minNodesNumAxis_0 - 1; i < max;)
	{
		for (size_t j = 0; j < minNodesNumAxis_0 - 1; j++)
		{
			_freeSimplexes.push_back(simp_buf = (new Simplex2())->GetPtrToUniquePtr());
			if (edge_buf = 
					FindEdge(
						**_freeNodes[i], 
						**_freeNodes[i + 1 + minNodesNumAxis_0]))
			{
				(*simp_buf)->edges[0] = edge_buf;
			}
			else
			{
				(*simp_buf)->edges[0] = 
					(new Edge2(
						**_freeNodes[i], 
						**_freeNodes[i + 1 + minNodesNumAxis_0]))
					->GetPtrToUniquePtr();

				_freeEdges.push_back((*simp_buf)->edges[0]);
			}
			if (edge_buf = 
					FindEdge(
						**_freeNodes[i + 1 + minNodesNumAxis_0], 
						**_freeNodes[i + minNodesNumAxis_0]))
			{
				(*simp_buf)->edges[1] = edge_buf;
			}
			else
			{
				(*simp_buf)->edges[1] = 
					(new Edge2(
						**_freeNodes[i + 1 + minNodesNumAxis_0], 
						**_freeNodes[i + minNodesNumAxis_0]))
					->GetPtrToUniquePtr();

				_freeEdges.push_back((*simp_buf)->edges[1]);
			}
			if (edge_buf = 
					FindEdge(
						**_freeNodes[i], 
						**_freeNodes[i + minNodesNumAxis_0]))
			{
				(*simp_buf)->edges[2] = edge_buf;
			}
			else
			{
				(*simp_buf)->edges[2] = 
					(new Edge2(
						**_freeNodes[i], 
						**_freeNodes[i + minNodesNumAxis_0]))
					->GetPtrToUniquePtr();

				_freeEdges.push_back((*simp_buf)->edges[2]);
			}

			(*(*simp_buf)->edges[0])->inclInSimplexes.push_back(simp_buf);
			(*(*simp_buf)->edges[1])->inclInSimplexes.push_back(simp_buf);
			(*(*simp_buf)->edges[2])->inclInSimplexes.push_back(simp_buf);

			for (auto &edge : (*simp_buf)->edges)
			{
				for (auto &node : (*edge)->nodes)
				{
					if (std::find(
						(*node)->inclInSimplexes.begin(),
						(*node)->inclInSimplexes.end(),
						(*simp_buf)->GetPtrToUniquePtr())
						== (*node)->inclInSimplexes.end())
					{
						(*node)->inclInSimplexes.push_back((*simp_buf)->GetPtrToUniquePtr());
					}
				}
			}

			i++;
		}
		i++;

		if (!(i < max))
		{
			break;
		}

		for (size_t j = 0; j < minNodesNumAxis_0; j++)
		{
			_freeSimplexes.push_back(simp_buf = (new Simplex2())->GetPtrToUniquePtr());
			if (edge_buf = 
					FindEdge(
						**_freeNodes[i], 
						**_freeNodes[i + 1 + minNodesNumAxis_0]))
			{
				(*simp_buf)->edges[0] = edge_buf;
			}
			else
			{
				(*simp_buf)->edges[0] = 
					(new Edge2(
						**_freeNodes[i], 
						**_freeNodes[i + 1 + minNodesNumAxis_0]))
					->GetPtrToUniquePtr();

				_freeEdges.push_back((*simp_buf)->edges[0]);
			}
			if (edge_buf = 
					FindEdge(
						**_freeNodes[i + 1 + minNodesNumAxis_0], 
						**_freeNodes[i + minNodesNumAxis_0]))
			{
				(*simp_buf)->edges[1] = edge_buf;
			}
			else
			{
				(*simp_buf)->edges[1] = 
					(new Edge2(
						**_freeNodes[i + 1 + minNodesNumAxis_0], 
						**_freeNodes[i + minNodesNumAxis_0]))
					->GetPtrToUniquePtr();

				_freeEdges.push_back((*simp_buf)->edges[1]);
			}
			if (edge_buf = FindEdge(**_freeNodes[i], **_freeNodes[i + minNodesNumAxis_0]))
			{
				(*simp_buf)->edges[2] = edge_buf;
			}
			else
			{
				(*simp_buf)->edges[2] = 
					(new Edge2(
						**_freeNodes[i], 
						**_freeNodes[i + minNodesNumAxis_0]))
					->GetPtrToUniquePtr();

				_freeEdges.push_back((*simp_buf)->edges[2]);
			}

			(*(*simp_buf)->edges[0])->inclInSimplexes.push_back(simp_buf);
			(*(*simp_buf)->edges[1])->inclInSimplexes.push_back(simp_buf);
			(*(*simp_buf)->edges[2])->inclInSimplexes.push_back(simp_buf);

			for (auto &edge : (*simp_buf)->edges)
			{
				for (auto &node : (*edge)->nodes)
				{
					if (std::find(
						(*node)->inclInSimplexes.begin(),
						(*node)->inclInSimplexes.end(),
						(*simp_buf)->GetPtrToUniquePtr())
						== (*node)->inclInSimplexes.end())
					{
						(*node)->inclInSimplexes.push_back((*simp_buf)->GetPtrToUniquePtr());
					}
				}
			}

			i++;
		}
		i++;
	}
}

void Polycrystalline2::GenerateUniformMesh()
{
	double polycrysSizeAxis[2];
	polycrysSizeAxis[0] = _maxShellNodesCoor[0] - _minShellNodesCoor[0];
	polycrysSizeAxis[1] = _maxShellNodesCoor[1] - _minShellNodesCoor[1];

	size_t minNodesNumAxis[2];
	minNodesNumAxis[0] = (size_t)(polycrysSizeAxis[0] / _l_av) + 3;
	minNodesNumAxis[1] = (size_t)(polycrysSizeAxis[1] / (0.5 * sqrt(3) * _l_av)) + 3;

	GenerateNodesEvenly(polycrysSizeAxis, minNodesNumAxis);
	GenerateSimplexesFromNodes(minNodesNumAxis[0]);
}

void Polycrystalline2::FitNodesToShellNodes()
{
	size_t maxi = _shellNodes.size();
	unique_ptr<Node2>* node_with_min_sqr_dist;
	double sqr_magn_buf;
	double min_sqr_dist;
	//#pragma omp parallel for private(node_with_min_sqr_dist, sqr_magn_buf, min_sqr_dist) firstprivate(maxi)
	for (size_t i = 0; i < maxi; i++)
	{
		node_with_min_sqr_dist = nullptr;
		min_sqr_dist = DBL_MAX;
		for (size_t j = 1, maxj = _freeNodes.size(); j < maxj; j++)
		{
			if (*_freeNodes[j])
			{
				sqr_magn_buf = (**_freeNodes[j] - **_shellNodes[i]).SqrMagnitude();
				if (sqr_magn_buf < min_sqr_dist)
				{
					min_sqr_dist = sqr_magn_buf;
					node_with_min_sqr_dist = _freeNodes[j];
				}
			}
		}

		//#pragma omp critical(FitFreeNodesToShellNodes)
		//{
		(*node_with_min_sqr_dist)->SetPosition((*_shellNodes[i])->GetPosition());
		(*node_with_min_sqr_dist)->belongsToShellNode = _shellNodes[i];
		for (auto &crys : (*_shellNodes[i])->inclInCryses)
		{
			if (std::find(
					(*node_with_min_sqr_dist)->belongsToCryses.begin(), 
					(*node_with_min_sqr_dist)->belongsToCryses.end(), 
					crys) 
				== (*node_with_min_sqr_dist)->belongsToCryses.end())
			{
				(*node_with_min_sqr_dist)->belongsToCryses.push_back(crys);
			}
		}
		//}
	}
}

void Polycrystalline2::FitNodesToShellEdges_OLD()
{
	size_t maxi = _freeEdges.size();
	//#pragma omp parallel for firstprivate(maxi)
	for (size_t i = 0; i < maxi; i++)
	{
		if (!(*_freeEdges[i]) ||
			(*(*_freeEdges[i])->nodes[0])->belongsToShellNode ||
			(*(*_freeEdges[i])->nodes[1])->belongsToShellNode ||
			(*(*_freeEdges[i])->nodes[0])->belongsToShellEdge ||
			(*(*_freeEdges[i])->nodes[1])->belongsToShellEdge)
		{
			continue;
		}

		for (size_t j = 0, maxj = _shellEdges.size(); j < maxj; j++)
		{
			Vector2 intersect_point;
			if (Vector2::SegmentsIntersection(
					intersect_point,
					(*(*_freeEdges[i])->nodes[0])->GetPosition(),
					(*(*_freeEdges[i])->nodes[1])->GetPosition(),
					(*(*_shellEdges[j])->nodes[0])->GetPosition(),
					(*(*_shellEdges[j])->nodes[1])->GetPosition()))
			{
				if (((*(*_freeEdges[i])->nodes[0])->GetPosition() - intersect_point).SqrMagnitude() <
					((*(*_freeEdges[i])->nodes[1])->GetPosition() - intersect_point).SqrMagnitude())
				{
					//#pragma omp critical(FitFreeNodesToShellEdges)
					//{
					(*(*_freeEdges[i])->nodes[0])->SetPosition(intersect_point);
					(*(*_freeEdges[i])->nodes[0])->belongsToShellEdge = _shellEdges[j];
					//}
				}
				else
				{
					//#pragma omp critical(FitFreeNodesToShellEdges)
					//{
					(*(*_freeEdges[i])->nodes[1])->SetPosition(intersect_point);
					(*(*_freeEdges[i])->nodes[1])->belongsToShellEdge = _shellEdges[j];
					//}
				}
			}
		}
	}
	
	//#pragma omp parallel for firstprivate(maxi)
	for (size_t i = 0; i < maxi; i++)
	{
		if (!(*_freeEdges[i]) ||
			!((*(*_freeEdges[i])->nodes[0])->belongsToShellEdge ||
			(*(*_freeEdges[i])->nodes[1])->belongsToShellEdge) ||
			(*(*_freeEdges[i])->nodes[0])->belongsToShellNode ||
			(*(*_freeEdges[i])->nodes[1])->belongsToShellNode ||
			((*(*_freeEdges[i])->nodes[0])->belongsToShellEdge &&
			 (*(*_freeEdges[i])->nodes[1])->belongsToShellEdge))
		{
			continue;
		}

		for (size_t j = 0, maxj = _shellEdges.size(); j < maxj; j++)
		{
			Vector2 intersect_point;
			if (Vector2::SegmentsIntersection(
					intersect_point,
					(*(*_freeEdges[i])->nodes[0])->GetPosition(),
					(*(*_freeEdges[i])->nodes[1])->GetPosition(),
					(*(*_shellEdges[j])->nodes[0])->GetPosition(),
					(*(*_shellEdges[j])->nodes[1])->GetPosition()))
			{
				if (!(*(*_freeEdges[i])->nodes[1])->belongsToShellEdge &&
					(*(*_freeEdges[i])->nodes[0])->belongsToShellEdge != _shellEdges[j])
				{
					//#pragma omp critical(FitFreeNodesToShellEdges)
					//{
					(*(*_freeEdges[i])->nodes[1])->SetPosition(intersect_point);
					(*(*_freeEdges[i])->nodes[1])->belongsToShellEdge = _shellEdges[j];
					//}
				}
				else if (!(*(*_freeEdges[i])->nodes[0])->belongsToShellEdge &&
					(*(*_freeEdges[i])->nodes[1])->belongsToShellEdge != _shellEdges[j])
				{
					//#pragma omp critical(FitFreeNodesToShellEdges)
					//{
					(*(*_freeEdges[i])->nodes[0])->SetPosition(intersect_point);
					(*(*_freeEdges[i])->nodes[0])->belongsToShellEdge = _shellEdges[j];
					//}
				}
			}
		}
	}

	size_t max = _freeNodes.size();
	//#pragma omp parallel for firstprivate(max)
	for (size_t i = 0; i < max; i++)
	{
		if ((*_freeNodes[i])->belongsToShellNode ||
			!(*_freeNodes[i])->belongsToShellEdge)
		{
			continue;
		}

		for (auto &crys : (*(*_freeNodes[i])->belongsToShellEdge)->inclInCryses)
		{
			if (std::find(
					(*_freeNodes[i])->belongsToCryses.begin(),
					(*_freeNodes[i])->belongsToCryses.end(),
					crys)
				== (*_freeNodes[i])->belongsToCryses.end())
			{
				(*_freeNodes[i])->belongsToCryses.push_back(crys);
			}
		}
	}
}

void Polycrystalline2::FitNodesToShellEdges()
{
	size_t nodes_num = _freeNodes.size();
	Vector2* shifts = new Vector2[nodes_num];

	int iters_num;
	std::cout << "Enter the iterations number: ";
	std::cin >> iters_num;
	for (int i = 0; i < iters_num; i++)
	{
		//#pragma omp parallel for firstprivate(nodes_num)
		for (size_t j = 0; j < nodes_num; j++)
		{
			if (*_freeNodes[j])
			{
				shifts[j] = ShiftToLavEdges(**_freeNodes[j]);
				shifts[j] += ShiftToFitMesh(**_freeNodes[j]);
			}
		}

		//#pragma omp parallel for firstprivate(nodes_num)
		for (size_t j = 0; j < nodes_num; j++)
		{
			if (*_freeNodes[j])
			{
				**_freeNodes[j] += shifts[j];
			}
		}
	}
}

void Polycrystalline2::FitMeshToShells()
{
	FitNodesToShellNodes();
	FitNodesToShellEdges();
}

void Polycrystalline2::DeleteExternalNodes()
{
	size_t maxi = _freeNodes.size();

	//#pragma omp parallel for firstprivate(maxi)
	for (size_t i = 0; i < maxi; i++)
	{
		if ((*_freeNodes[i])->belongsToShellEdge ||
			(*_freeNodes[i])->belongsToShellNode)
		{
			continue;
		}

		for (auto &crys : crystallites)
		{
			if (crys->Contains(**_freeNodes[i]))
			{
				(*_freeNodes[i])->belongsToCryses.push_back(crys);
			}
		}
	}

	for (auto &node : _freeNodes)
	{
		if (*node && (*node)->belongsToCryses.empty())
		{
			delete node->release();
		}
	}
}

void Polycrystalline2::DeleteFarNodes()
{
	size_t nodes_num = _freeNodes.size();
	bool* is_far_nodes = new bool[nodes_num];
	bool* is_inside_nodes = new bool[nodes_num];
	
	//#pragma omp parallel for firstprivate(nodes_num)
	for (size_t i = 0; i < nodes_num; i++)
	{
		is_inside_nodes[i] = false;

		if (!*_freeNodes[i] ||
			(*_freeNodes[i])->belongsToShellNode)
		{
			is_inside_nodes[i] = true;
			continue;
		}

		for (auto &crys : crystallites)
		{
			if (crys->Contains(**_freeNodes[i]))
			{
				is_inside_nodes[i] = true;
				break;
			}
		}
	}

	//#pragma omp parallel for firstprivate(nodes_num)
	for (size_t i = 0; i < nodes_num; i++)
	{
		is_far_nodes[i] = false;

		if (*_freeNodes[i])
		{
			bool do_continue = false;

			double min_sqr_dist_to_s_node = DBL_MAX;
			double sqr_magn_buf;
			for (auto &s_node : _shellNodes)
			{
				if ((**s_node - **_freeNodes[i]).SqrMagnitude() < 4.0 * _l_av * _l_av)
				{
					is_far_nodes[i] = false;
					do_continue = true;
					break;
				}
			}
			if (do_continue)
			{
				continue;
			}
			

			Vector2 proj;
			for (auto &s_edge : _shellEdges)
			{
				if (*s_edge &&
					Vector2::Project(
						proj,
						(*_freeNodes[i])->GetPosition(),
						(*(*s_edge)->nodes[0])->GetPosition(),
						(*(*s_edge)->nodes[1])->GetPosition()))
				{
					if ((proj - (*_freeNodes[i])->GetPosition()).SqrMagnitude() < 4.0 * _l_av * _l_av)
					{
						is_far_nodes[i] = false;
						do_continue = true;
						break;
					}
				}
			}
			if (do_continue)
			{
				continue;
			}

			is_far_nodes[i] = true;
		}
	}

	for (size_t i = 0; i < nodes_num; i++)
	{
		if (!is_inside_nodes[i] &&
			is_far_nodes[i])
		{
			delete _freeNodes[i]->release();
		}
	}

	delete is_far_nodes;
	delete is_inside_nodes;
}

void Polycrystalline2::DivideExtendedSimplexes()
{
	//#pragma omp parallel for firstprivate(simps_num)
	for (auto &simp : _freeSimplexes)
	{
		if (*simp)
		{
			(*simp)->UpdateAverageEdgesLength();
		}
	}

	for (size_t i = 0, max = _freeEdges.size(); i < max; i++)
	{
		if (*_freeEdges[i])
		{
			for (auto &simp : (*_freeEdges[i])->inclInSimplexes)
			{
				if ((*_freeEdges[i])->SqrMagnitude() > 1.937 * (*simp)->averageEdgesLength * (*simp)->averageEdgesLength)
				{
					(*_freeEdges[i])->MakeTwoInstead(_freeSimplexes, _freeEdges, _freeNodes);
					break;
				}
				if ((*(*simp)->MinEdge())->SqrMagnitude() < 0.375 * 0.375 * (*simp)->averageEdgesLength * (*simp)->averageEdgesLength)
				{
					(*(*simp)->MaxEdge())->MakeTwoInstead(_freeSimplexes, _freeEdges, _freeNodes);
					break;
				}
			}
		}
	}
}

void Polycrystalline2::DivideCrossingEdges()
{
	for (size_t i = 0, max = _freeEdges.size(); i < max; i++)
	{
		if (*_freeEdges[i] &&
			(*_freeEdges[i])->SqrMagnitude() > 0.25 * _l_av * _l_av &&
			(*_freeEdges[i])->IntersectsWith(_shellEdges))
		{
			(*_freeEdges[i])->MakeTwoInstead(_freeSimplexes, _freeEdges, _freeNodes);
		}
	}
}

void Polycrystalline2::MinMaxEdges(double& min, double& max)
{
	size_t edges_num = _freeEdges.size();

	min = ULLONG_MAX;
	max = 0;
	double sqr_length_buf = -1.0;
	for (size_t i = 0; i < edges_num; i++)
	{
		if (*_freeEdges[i] &&
			(*_freeEdges[i])->nodes[0] &&
			(*_freeEdges[i])->nodes[1])
		{
			sqr_length_buf = (**(*_freeEdges[i])->nodes[0] - **(*_freeEdges[i])->nodes[1]).SqrMagnitude();
			if (sqr_length_buf > max)
			{
				max = sqr_length_buf;
			}
			if (sqr_length_buf < min)
			{
				min = sqr_length_buf;
			}
		}
	}

	min = sqrt(min);
	max = sqrt(max);
}

Vector2 Polycrystalline2::ShiftToFitMesh(const Node2& node)
{
	Vector2 shift;
	Vector2 proj;
	const double c = 0.2;
	for (auto &s_edge : _shellEdges)
	{
		if (*s_edge &&
			Vector2::Project(
				proj,
				node.GetPosition(),
				(*(*s_edge)->nodes[0])->GetPosition(),
				(*(*s_edge)->nodes[1])->GetPosition()))
		{
			Vector2 to_proj = proj - node.GetPosition();
			Vector2 buf;
			if (to_proj.SqrMagnitude() > _l_av * _l_av * 0.25)
			{
				buf = Vector2(0.0, 0.0);
			}
			else if (to_proj.SqrMagnitude() < _l_av * 0.2)
			{
				buf = to_proj;
			}
			else
			{
				double sqr_dist = to_proj.SqrMagnitude();
				buf = c * to_proj * _l_av / (sqr_dist);
				if (buf.SqrMagnitude() > _l_av * _l_av * 0.0025)
				{
					buf = buf.Normalize() * _l_av * 0.05;
				}
			}

			shift += buf;
		}
	}

	return shift;
}

double Sign(double val)
{
	return val > 0.0 ? 1.0 : -1.0;
}

Vector2 Polycrystalline2::ShiftToLavEdges(const Node2& node)
{
	Vector2 shift;
	const double b = 0.15;
	for (auto &neighbor : node.neighbors)
	{
		if (*neighbor)
		{
			Vector2& vec = **neighbor - node;
			double buf = b * (1.0 - _l_av / vec.Magnitude());
			if ((*neighbor)->belongsToShellNode)
			{
				shift += vec * 2 * (buf);
			}
			else if ((*neighbor)->belongsToShellEdge)
			{
				double k = 2.0 - Vector2::Cos(vec, **(*(*neighbor)->belongsToShellEdge)->nodes[0] - **(*(*neighbor)->belongsToShellEdge)->nodes[1]);
				shift += vec * k * (buf);
			}
			else
			{
				shift += vec * (buf);
			}
		}
	}

	return shift;
}

Vector2 Polycrystalline2::ShiftToEquilateralSimplexes(const Node2& node)
{
	Vector2 shift;
	const double a = 0.002;
	for (auto &simp : node.inclInSimplexes)
	{
		if (*simp)
		{
			for (auto &neighbor : node.neighbors)
			{
				if ((*simp)->IsContaining(**neighbor))
				{
					Vector2& vec = **neighbor - node;
					double buf = 1.0 - (*simp)->averageEdgesLength / vec.Magnitude();
					if ((*neighbor)->belongsToShellNode)
					{
						shift += vec * 2.0 * a * (buf + 2.0 * buf * buf * Sign(buf));
					}
					else if ((*neighbor)->belongsToShellEdge)
					{
						double k = 2.0 - Vector2::Cos(vec, **(*(*neighbor)->belongsToShellEdge)->nodes[0] - **(*(*neighbor)->belongsToShellEdge)->nodes[1]);
						shift += vec * k * a * (buf + 2.0 * buf * buf * Sign(buf));
					}
					else
					{
						shift += vec * a * (buf + 2.0 * buf * buf * Sign(buf));
					}
				}
			}
		}
	}

	return shift;
}

void Polycrystalline2::DistributeNodesEvenly()
{
	size_t nodes_num = _freeNodes.size();
	Vector2* shifts = new Vector2[nodes_num];

	//int continue_flag = 1;
	//double sufficient_delta = (_l_max - _l_min) * 0.001;
	//double prev_min, prev_max;
	//double min, max;
	//double d_min, d_max;
	int iters_num;
	std::cout << "Enter the iterations number: ";
	std::cin >> iters_num;
	//MinMaxEdges(min, max);
	for (int i = 0; i < iters_num; i++)
	{
		//prev_min = min;
		//prev_max = max;

		for (auto &simp : _freeSimplexes)
		{
			if (*simp)
			{
				(*simp)->UpdateAverageEdgesLength();
			}
		}

		//#pragma omp parallel for firstprivate(nodes_num)
		for (size_t j = 0; j < nodes_num; j++)
		{
			if (*_freeNodes[j])
			{
				shifts[j] = ShiftToLavEdges(**_freeNodes[j]);
			}
		}
		//#pragma omp parallel for firstprivate(nodes_num)
		for (size_t j = 0; j < nodes_num; j++)
		{
			if (*_freeNodes[j])
			{
				**_freeNodes[j] += shifts[j];
			}
		}

		//MinMaxEdges(min, max);
		//d_min = min - prev_min;
		//d_max = max - prev_max;
		
		//if (continue_flag == 1 &&
		//	abs(d_min) < sufficient_delta &&
		//	abs(d_max) < sufficient_delta)
		//{
		//	std::cout << "Stoped after " << iter << " iterations.\n"
		//				 "Average delta at the end: " << (abs(d_min) + abs(d_max)) * 0.5;
		//	std::cout << "\nEnter '0' to continue and '1' to stop\n";
		//	std::cin >> continue_flag;
		//	if (continue_flag == 1)
		//	{
		//		continue_flag = 0;
		//		break;
		//	}
		//}
	}
	//if (continue_flag == 1)
	//{
	//	std::cout << "Stoped after " << iterations << " iterations.\n"
	//		"Average delta at the end: " << (abs(d_min) + abs(d_max)) * 0.5 << '\n';
	//}

	delete[] shifts;
}

void Polycrystalline2::MakeSimplexesEquilateral()
{
	size_t nodes_num = _freeNodes.size();
	Vector2* shifts = new Vector2[nodes_num];

	int iters_num;
	std::cout << "Enter the iterations number: ";
	std::cin >> iters_num;
	for (int i = 0; i < iters_num; i++)
	{
		for (auto &simp : _freeSimplexes)
		{
			if (*simp)
			{
				(*simp)->UpdateAverageEdgesLength();
			}
		}

		//#pragma omp parallel for firstprivate(nodes_num)
		for (size_t i = 0; i < nodes_num; i++)
		{
			if (*_freeNodes[i])
			{
				shifts[i] = ShiftToEquilateralSimplexes(**_freeNodes[i]);
			}
		}
		//#pragma omp parallel for firstprivate(nodes_num)
		for (size_t i = 0; i < nodes_num; i++)
		{
			if (*_freeNodes[i])
			{
				**_freeNodes[i] += shifts[i];
			}
		}
	}

	delete[] shifts;
}

const bool Polycrystalline2::Contains(const Node2& node) const
{
	for (auto &crys : crystallites)
	{
		if (crys->Contains(node))
		{
			return true;
		}
	}

	return false;
}

void AddNodesData(ofstream& nodesData, const vector<unique_ptr<Node2>*>& nodes)
{
	bool is_first = true;
	size_t index = 0;
	for (size_t i = 0, max = nodes.size(); i < max; i++)
	{
		if (*nodes[i])
		{
			if (!is_first)
			{
				nodesData << '\n';
			}
			nodesData << (*nodes[i])->GetPosition()[0] << ' ';
			nodesData << (*nodes[i])->GetPosition()[1];
			(*nodes[i])->globalNum = index++;
			is_first = false;
		}
	}
}

void AddNodesData(vector<double>& nodesData, const vector<unique_ptr<Node2>*>& nodes)
{
	size_t index = 0;
	for (size_t i = 0, max = nodes.size(); i < max; i++)
	{
		if (*nodes[i])
		{
			nodesData.push_back((*nodes[i])->GetPosition()[0]);
			nodesData.push_back((*nodes[i])->GetPosition()[1]);
			(*nodes[i])->globalNum = index++;
		}
	}
}

// Different for 3 dimensions.
void AddFENodesDataFromSimpex(ofstream& fenData, const Simplex2& simp, bool isFirst)
{
	if (!isFirst)
	{
		fenData << '\n';
	}
	fenData << (*(*simp.edges[0])->nodes[0])->globalNum << ' ';
	fenData << (*(*simp.edges[0])->nodes[1])->globalNum << ' ';
	if ((*simp.edges[1])->nodes[0] == (*simp.edges[0])->nodes[0] ||
		(*simp.edges[1])->nodes[0] == (*simp.edges[0])->nodes[1])
	{
		fenData << (*(*simp.edges[1])->nodes[1])->globalNum;
	}
	else
	{
		fenData << (*(*simp.edges[1])->nodes[0])->globalNum;
	}
}

// Different for 3 dimensions.
void AddFENodesDataFromSimpex(vector<size_t>& fenData, const Simplex2& simp)
{
	fenData.push_back((*(*simp.edges[0])->nodes[0])->globalNum);
	fenData.push_back((*(*simp.edges[0])->nodes[1])->globalNum);
	if ((*simp.edges[1])->nodes[0] == (*simp.edges[0])->nodes[0] ||
		(*simp.edges[1])->nodes[0] == (*simp.edges[0])->nodes[1])
	{
		fenData.push_back((*(*simp.edges[1])->nodes[1])->globalNum);
	}
	else
	{
		fenData.push_back((*(*simp.edges[1])->nodes[0])->globalNum);
	}
}

void AddFENodesData(ofstream& feNodesData, const list<unique_ptr<Simplex2>*>& simplexes)
{
	bool is_first = true;
	for (list<unique_ptr<Simplex2>*>::const_iterator simp_iter = simplexes.begin();
		simp_iter != simplexes.end();
		simp_iter++)
	{
		if (**simp_iter)
		{
			AddFENodesDataFromSimpex(feNodesData, ***simp_iter, is_first);
			is_first = false;
		}
	}
}

void AddFENodesData(vector<size_t>& feNodesData, const list<unique_ptr<Simplex2>*>& simplexes)
{
	for (list<unique_ptr<Simplex2>*>::const_iterator simp_iter = simplexes.begin();
		simp_iter != simplexes.end();
		simp_iter++)
	{
		if (**simp_iter)
		{
			AddFENodesDataFromSimpex(feNodesData, ***simp_iter);
		}
	}
}

// Different for 3 dimensions.
void Polycrystalline2::InputData(ifstream& shell_nodes, ifstream& cryses_edges, ifstream& gener_params)
{
	_minShellNodesCoor[0] = _minShellNodesCoor[1] = DBL_MAX;
	_maxShellNodesCoor[0] = _maxShellNodesCoor[1] = DBL_MIN;

	double dbuf[2];
	while (!shell_nodes.eof())
	{
		shell_nodes >> dbuf[0];
		if (dbuf[0] < _minShellNodesCoor[0])
		{
			_minShellNodesCoor[0] = dbuf[0];
		}
		if (dbuf[0] > _maxShellNodesCoor[0])
		{
			_maxShellNodesCoor[0] = dbuf[0];
		}
		shell_nodes >> dbuf[1];
		if (dbuf[1] < _minShellNodesCoor[1])
		{
			_minShellNodesCoor[1] = dbuf[1];
		}
		if (dbuf[1] > _maxShellNodesCoor[1])
		{
			_maxShellNodesCoor[1] = dbuf[1];
		}
		_shellNodes.push_back(
			(new ShellNode2(
				dbuf[0], 
				dbuf[1]))
			->GetPtrToUniquePtr());
	}
	shell_nodes.clear();
	shell_nodes.seekg(0);

	size_t ibuf[2];
	size_t nums_num;
	unique_ptr<ShellEdge2>* shell_edge_buf = nullptr;
	for (list<Crystallite2*>::iterator crys_iter; !cryses_edges.eof();)
	{
		if (crystallites.empty())
		{
			crystallites.push_back(new Crystallite2());
			crys_iter = crystallites.begin();
		}
		else
		{
			crystallites.push_back(new Crystallite2());
			crys_iter++;
		}
		cryses_edges >> nums_num;
		for (size_t j = 0; j < nums_num; j++)
		{
			cryses_edges >> ibuf[0];
			cryses_edges >> ibuf[1];
			
			for (auto &s_edge : (*_shellNodes[ibuf[0]])->inclInEdges)
			{
				if ((*s_edge)->IsContaining(**_shellNodes[ibuf[1]]))
				{
					shell_edge_buf = s_edge;
				}
			}

			if (shell_edge_buf)
			{
				(*crys_iter)->shellEdges.push_back(shell_edge_buf);
				(*shell_edge_buf)->inclInCryses.push_back(*crys_iter);
			}
			else
			{
				(*crys_iter)->shellEdges.push_back(
					shell_edge_buf =
						(new ShellEdge2(
							**_shellNodes[ibuf[0]],
							**_shellNodes[ibuf[1]]))->GetPtrToUniquePtr());

				(*_shellNodes[ibuf[0]])->inclInCryses.push_back(*crys_iter);
				(*_shellNodes[ibuf[1]])->inclInCryses.push_back(*crys_iter);
				(*shell_edge_buf)->inclInCryses.push_back(*crys_iter);
				_shellEdges.push_back(shell_edge_buf);
			}
			shell_edge_buf = nullptr;
		}
	}
	cryses_edges.clear();
	cryses_edges.seekg(0);

	gener_params >> _l_min;
	gener_params >> _l_max;
	gener_params.clear();
	gener_params.seekg(0);

	_l_av = (_l_min + _l_max) * 0.5;
}

// Different for 3 dimensions.
void Polycrystalline2::InputData(const vector<double>& shell_nodes, const vector<size_t>& cryses_edges, const vector<size_t>& gener_params)
{
	_minShellNodesCoor[0] = _minShellNodesCoor[1] = DBL_MAX;
	_maxShellNodesCoor[0] = _maxShellNodesCoor[1] = DBL_MIN;

	for (size_t i = 0, max = shell_nodes.size(); i < max; i += 2) // i += 3 for 3 dimensions
	{
		if (shell_nodes[i] < _minShellNodesCoor[0])
		{
			_minShellNodesCoor[0] = shell_nodes[i];
		}
		if (shell_nodes[i] > _maxShellNodesCoor[0])
		{
			_maxShellNodesCoor[0] = shell_nodes[i];
		}
		if (shell_nodes[i + 1] < _minShellNodesCoor[1])
		{
			_minShellNodesCoor[1] = shell_nodes[i + 1];
		}
		if (shell_nodes[i + 1] > _maxShellNodesCoor[1])
		{
			_maxShellNodesCoor[1] = shell_nodes[i + 1];
		}

		_shellNodes.push_back(
			(new ShellNode2(
				shell_nodes[i],
				shell_nodes[i + 1]))
			->GetPtrToUniquePtr());

		size_t ibuf[2];
		size_t nums_num;
		unique_ptr<ShellEdge2>* shell_edge_buf;
		vector<size_t>::const_iterator cp_iter = cryses_edges.begin();
		for (list<Crystallite2*>::iterator crys_iter; cp_iter != cryses_edges.end(); crys_iter++)
		{
			if (crystallites.empty())
			{
				crystallites.push_back(new Crystallite2());
				crys_iter = crystallites.begin();
			}
			else
			{
				crystallites.push_back(new Crystallite2());
				crys_iter++;
			}
			nums_num = *(cp_iter++);
			for (size_t j = 0; j < nums_num; j++)
			{
				for (auto &s_edge : (*_shellNodes[ibuf[0]])->inclInEdges)
				{
					if ((*s_edge)->IsContaining(**_shellNodes[ibuf[1]]))
					{
						shell_edge_buf = s_edge;
					}
				}
				if (shell_edge_buf)
				{
					(*crys_iter)->shellEdges.push_back(shell_edge_buf);
					(*shell_edge_buf)->inclInCryses.push_back(*crys_iter);
				}
				else
				{
					(*crys_iter)->shellEdges.push_back(
						shell_edge_buf =
							(new ShellEdge2(
								**_shellNodes[*(cp_iter)],
								**_shellNodes[*(cp_iter + 1)]))
						->GetPtrToUniquePtr());

					(*_shellNodes[*(cp_iter)])->inclInCryses.push_back(*crys_iter);
					(*_shellNodes[*(cp_iter + 1)])->inclInCryses.push_back(*crys_iter);
					(*shell_edge_buf)->inclInCryses.push_back(*crys_iter);
					_shellEdges.push_back(shell_edge_buf);
				}
				shell_edge_buf = nullptr;

				cp_iter += 2;
			}
		}
	}
	_l_min = gener_params[0];
	_l_max = gener_params[1];
	_l_av = (_l_min + _l_max) * 0.5;
}

void Polycrystalline2::OutputData(ofstream& nodesData, ofstream& feNodesData) const
{
	AddNodesData(nodesData, _freeNodes);
	AddFENodesData(feNodesData, _freeSimplexes);
}

void Polycrystalline2::OutputData(vector<double>& nodesData, vector<size_t>& feNodesData) const
{
	AddNodesData(nodesData, _freeNodes);
	AddFENodesData(feNodesData, _freeSimplexes);
}

Polycrystalline2::Polycrystalline2() {}

Polycrystalline2::Polycrystalline2(ifstream& shell_nodes, ifstream& cryses_edges, ifstream& gener_params)
{
	InputData(shell_nodes, cryses_edges, gener_params);
}

Polycrystalline2::Polycrystalline2(const vector<double>& shell_nodes, const vector<size_t>& cryses_edges, const vector<size_t>& gener_params)
{
	InputData(shell_nodes, cryses_edges, gener_params);
}

Polycrystalline2::~Polycrystalline2()
{
	for (auto &crys : crystallites)
	{
		if (crys)
		{
			delete crys;
			crys = nullptr;
		}
	}

	for (auto &simp : _freeSimplexes)
	{
		if (*simp)
		{
			delete simp->release();
		}
	}

	for (auto &ptr : _shellEdges)
	{
		delete ptr;
	}
	for (auto &ptr : _shellNodes)
	{
		delete ptr;
	}
	for (auto &ptr : _freeSimplexes)
	{
		delete ptr;
	}
	for (auto &ptr : _freeEdges)
	{
		delete ptr;
	}
	for (auto &ptr : _freeNodes)
	{
		delete ptr;
	}
}