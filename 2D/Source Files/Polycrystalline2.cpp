#include "Polycrystalline2.h"
#include <algorithm>
#include <iostream>


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

	throw std::exception("Error in function Polycrystalline2::FindEdge");
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
		for (size_t j = 0, maxj = _freeNodes.size(); j < maxj; j++)
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
		(*node_with_min_sqr_dist)->belongsToCrys = true;
		//}
	}
}

void Polycrystalline2::AttachNodesToShellEdges(size_t iterations_num)
{
	size_t shell_edges_num = _shellEdges.size();

	//#pragma omp parallel firstprivate(shell_edges_num)
	//{
	//#pragma omp for
	for (size_t i = 0; i < shell_edges_num; i++)
	{
		(*_shellEdges[i])->AttachNodes(_freeNodes);
	}
	//#pragma omp for
	for (size_t i = 0; i < shell_edges_num; i++)
	{
		if ((*_shellEdges[i])->attachedNodes.size() == 0)
		{
			continue;
		}

		size_t j = 0;
		for (size_t k = 0; k < shell_edges_num; k++)
		{
			if (k == i)
			{
				continue;
			}

			if ((*_shellEdges[k])->ContainsAttachedNode((*_shellEdges[i])->attachedNodes[j]))
			{
				(*_shellEdges[i])->ChangeAttachedNode(j, _freeNodes);
			}
		}

		j = (*_shellEdges[i])->attachedNodes.size() - 1;
		if (j == 0)
		{
			continue;
		}

		for (size_t k = 0; k < shell_edges_num; k++)
		{
			if (k == i)
			{
				continue;
			}

			if ((*_shellEdges[k])->ContainsAttachedNode((*_shellEdges[i])->attachedNodes[j]))
			{
				(*_shellEdges[i])->ChangeAttachedNode(j, _freeNodes);
			}
		}
	}
	//#pragma omp for
	for (size_t i = 0; i < shell_edges_num; i++)
	{
		(*_shellEdges[i])->SetAttachedNodesStartVectorsToEdge();
	}
	//}

	DeleteExternalNodes();

	for (size_t i = 1; i < iterations_num + 1; i++)
	{
		//#pragma omp parallel for firstprivate(shell_edges_num)
		for (size_t j = 0; j < shell_edges_num; j++)
		{
			(*_shellEdges[j])->SetAttachedNodesDistanceFromStartPositionToEdge(i, iterations_num);
		}
		DistributeNodesEvenly(1);
	}
}

void Polycrystalline2::FitNodesToShellEdges()
{
	size_t nodes_num = _freeNodes.size();
	size_t shell_edges_num = _shellEdges.size();
	Vector2* shifts = new Vector2[nodes_num];

	int iters_num;
	std::cout << "Enter the iterations number: ";
	std::cin >> iters_num;
	for (int i = 0; i < iters_num; i++)
	{
		//#pragma omp parallel firstprivate(nodes_num)
		//{
		//#pragma omp for
		for (size_t j = 0; j < nodes_num; j++)
		{
			if (*_freeNodes[j])
			{
				shifts[j] = ShiftToLavEdges(**_freeNodes[j]);
				shifts[j] += ShiftToFitMesh(**_freeNodes[j]);
			}
		}
		//#pragma omp for
		for (size_t j = 0; j < nodes_num; j++)
		{
			if (*_freeNodes[j])
			{
				**_freeNodes[j] += shifts[j];
			}
		}
		//}
	}

	AttachNodesToShellEdges(iters_num);
}

void Polycrystalline2::FitMeshToShells()
{
	FitNodesToShellNodes();
	FitNodesToShellEdges();
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

void Polycrystalline2::ErasePtrsToNullptrFromVectors()
{
	//#pragma omp parallel num_threads(3)
	//{
	//#pragma omp single
	ErasePtrsToNullptr(_freeNodes);
	//#pragma omp single
	ErasePtrsToNullptr(_freeEdges);
	//#pragma omp single
	_freeSimplexes.remove_if([](unique_ptr<Simplex2>*& ptr) { return *ptr ? false : true; });
	//}
}

void Polycrystalline2::DeleteExternalNodes()
{
	size_t nodes_num = _freeNodes.size();

	//#pragma omp parallel for firstprivate(nodes_num)
	for (size_t i = 0; i < nodes_num; i++)
	{
		if ((*_freeNodes[i])->belongsToShellEdge ||
			(*_freeNodes[i])->belongsToShellNode)
		{
			(*_freeNodes[i])->belongsToCrys = true;
			continue;
		}

		for (auto &crys : crystallites) // Check only for external polycrystalline shell
		{
			if (crys->Contains(**_freeNodes[i]))
			{
				(*_freeNodes[i])->belongsToCrys = true;
				break;
			}
		}
	}

	for (auto &node : _freeNodes)
	{
		if (*node && !(*node)->belongsToCrys)
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
	
	//#pragma omp parallel firstprivate(nodes_num)
	//{
	//#pragma omp for
	for (size_t i = 0; i < nodes_num; i++)
	{
		is_inside_nodes[i] = false;

		if (!*_freeNodes[i] ||
			(*_freeNodes[i])->belongsToShellNode)
		{
			is_inside_nodes[i] = true;
			continue;
		}

		for (auto &crys : crystallites) // Check only for external polycrystalline shell
		{
			if (crys->Contains(**_freeNodes[i]))
			{
				is_inside_nodes[i] = true;
				break;
			}
		}
	}

	//#pragma omp for
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
	//}

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
				buf[0] = 0.0;
				buf[1] = 0.0;
			}
			else
			{
				double sqr_dist = to_proj.SqrMagnitude();
				buf = c * to_proj * _l_av / (sqr_dist);
			}

			shift += buf;
		}
	}

	if (shift.SqrMagnitude() > _l_av * _l_av * 0.0025)
	{
		shift = shift.Normalize() * _l_av * 0.05;
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
			Vector2 vec = **neighbor - node;
			double buf = b * (1.0 - _l_av / vec.Magnitude());
			if ((*neighbor)->belongsToShellNode)
			{
				shift += vec * 2.0 * buf;
			}
			else if ((*neighbor)->belongsToShellEdge)
			{
				double k = 2.0 - abs(Vector2::Cos(vec, **(*(*neighbor)->belongsToShellEdge)->nodes[0] - **(*(*neighbor)->belongsToShellEdge)->nodes[1]));
				shift += vec * k * buf;
			}
			else
			{
				shift += vec * buf;
			}
		}
	}

	if (shift.SqrMagnitude() > _l_av * _l_av * 0.0025)
	{
		shift = shift.Normalize() * _l_av * 0.05;
	}

	return shift;
}

//Vector2 Polycrystalline2::ShiftToEquilateralSimplexes(const Node2& node)
//{
//	Vector2 shift;
//	const double a = 0.002;
//	for (auto &simp : node.inclInSimplexes)
//	{
//		if (*simp)
//		{
//			for (auto &neighbor : node.neighbors)
//			{
//				if ((*simp)->IsContaining(**neighbor))
//				{
//					Vector2& vec = **neighbor - node;
//					double buf = 1.0 - (*simp)->averageEdgesLength / vec.Magnitude();
//					if ((*neighbor)->belongsToShellNode)
//					{
//						shift += vec * 2.0 * a * (buf + 2.0 * buf * buf * Sign(buf));
//					}
//					else if ((*neighbor)->belongsToShellEdge)
//					{
//						double k = 2.0 - Vector2::Cos(vec, **(*(*neighbor)->belongsToShellEdge)->nodes[0] - **(*(*neighbor)->belongsToShellEdge)->nodes[1]);
//						shift += vec * k * a * (buf + 2.0 * buf * buf * Sign(buf));
//					}
//					else
//					{
//						shift += vec * a * (buf + 2.0 * buf * buf * Sign(buf));
//					}
//				}
//			}
//		}
//	}
//
//	return shift;
//}

void Polycrystalline2::DistributeNodesEvenly(size_t iterations_num)
{
	size_t nodes_num = _freeNodes.size();
	Vector2* shifts = new Vector2[nodes_num];

	if (iterations_num == 0)
	{
		std::cout << "Enter the iterations number: ";
		std::cin >> iterations_num;
	}
	for (int i = 0; i < iterations_num; i++)
	{
		//for (auto &simp : _freeSimplexes)
		//{
		//	if (*simp)
		//	{
		//		(*simp)->UpdateAverageEdgesLength();
		//	}
		//}

		//#pragma omp parallel firstprivate(nodes_num)
		//{
		//#pragma omp for
		for (size_t j = 0; j < nodes_num; j++)
		{
			if (*_freeNodes[j])
			{
				shifts[j] = ShiftToLavEdges(**_freeNodes[j]);
			}
		}
		//#pragma omp for
		for (size_t j = 0; j < nodes_num; j++)
		{
			if (*_freeNodes[j])
			{
				**_freeNodes[j] += shifts[j];
			}
		}
		//}
	}

	delete[] shifts;
}

void Polycrystalline2::DelaunayPostprocessing()
{
	unique_ptr<Node2>* around_nodes[2];
	size_t edges_num = _freeEdges.size();
	for (size_t i = 0; i < edges_num; i++)
	{
		if (!*_freeEdges[i] ||
			(*(*_freeEdges[i])->nodes[0])->belongsToShellEdge == (*(*_freeEdges[i])->nodes[1])->belongsToShellEdge ||
			(*(*_freeEdges[i])->nodes[0])->belongsToShellEdge && (*(*_freeEdges[i])->nodes[1])->belongsToShellNode ||
			(*(*_freeEdges[i])->nodes[1])->belongsToShellEdge && (*(*_freeEdges[i])->nodes[0])->belongsToShellNode ||
			(*_freeEdges[i])->inclInSimplexes.size() < 2)
		{
			continue;
		}
		if ((*(*_freeEdges[i])->nodes[0])->belongsToShellNode)
		{
			bool continue_ = false;
			for (auto &s_edge : (*(*(*_freeEdges[i])->nodes[0])->belongsToShellNode)->inclInEdges)
			{
				if (*s_edge &&
					(*s_edge)->ContainsAttachedNode((*_freeEdges[i])->nodes[0]))
				{
					continue_ = true;
					break;
				}
			}
			if (continue_)
			{
				continue;
			}
		}
		
		around_nodes[0] = (*(*_freeEdges[i])->inclInSimplexes.front())->FindNodeNotIncludedInEdge(**_freeEdges[i]);
		around_nodes[1] = (*(*_freeEdges[i])->inclInSimplexes.back())->FindNodeNotIncludedInEdge(**_freeEdges[i]);
		if ((*around_nodes[0])->belongsToShellEdge && (*around_nodes[1])->belongsToShellNode ||
			(*around_nodes[1])->belongsToShellEdge && (*around_nodes[0])->belongsToShellNode)
		{
			continue;
		}

		(*_freeEdges[i])->FlipIfNeeded(_freeEdges, _freeSimplexes);
	}
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
				(*shell_edge_buf)->inclInCrysesNum++;
			}
			else
			{
				(*crys_iter)->shellEdges.push_back(
					shell_edge_buf =
						(new ShellEdge2(
							**_shellNodes[ibuf[0]],
							**_shellNodes[ibuf[1]]))->GetPtrToUniquePtr());

				(*shell_edge_buf)->inclInCrysesNum++;
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
					(*shell_edge_buf)->inclInCrysesNum++;
				}
				else
				{
					(*crys_iter)->shellEdges.push_back(
						shell_edge_buf =
							(new ShellEdge2(
								**_shellNodes[*(cp_iter)],
								**_shellNodes[*(cp_iter + 1)]))
						->GetPtrToUniquePtr());

					(*shell_edge_buf)->inclInCrysesNum++;
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