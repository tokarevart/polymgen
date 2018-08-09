#include "Polycrystalline2.h"
#include <algorithm>

#define DET(a, b, c, d) \
		(a * d - b * c)
#define BETWEEN(p0_coor, p1_coor, p) \
		(std::min(p0_coor, p1_coor) - DBL_EPSILON < p && p < std::max(p0_coor, p1_coor) + DBL_EPSILON)

void Polycrystalline2::Debug()
{
	//Crystallite2* crys_buf;
	//crystallites.push_back(crys_buf = new Crystallite2());
	(*crystallites.begin())->simplexes = std::move(_freeSimplexes);
}

void Polycrystalline2::GenerateFreeNodesEvenly(double* const polycrysSizeAxis, size_t* const minNodesNumAxis)
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
			new Node2(
				node_coors[0],
				node_coors[1]));
		i++;
		for (size_t j = 0; j < minNodesNumAxis[0]; j++)
		{
			node_coors[0] += _l_av;
			_freeNodes.push_back(
				new Node2(
					node_coors[0],
					node_coors[1]));
			i++;
		}

		if (!(i < freeNodesNum))
		{
			break;
		}

		node_coors[0] -= minNodesNumAxis[0] * _l_av - _l_av * 0.5;
		node_coors[1] += delta_node_coors_1;
		_freeNodes.push_back(
			new Node2(
				node_coors[0],
				node_coors[1]));
		i++;
		for (size_t j = 0; j < minNodesNumAxis[0] - 1; j++)
		{
			node_coors[0] += _l_av;
			_freeNodes.push_back(
				new Node2(
					node_coors[0],
					node_coors[1]));
			i++;
		}
		node_coors[0] -= minNodesNumAxis[0] * _l_av - _l_av * 0.5;
		node_coors[1] += delta_node_coors_1;
	}
}

unique_ptr<Edge2>* Polycrystalline2::FindFreeEdge(const Node2& node0, const Node2& node1)
{
	if (std::find(node0.neighbors.begin(), node0.neighbors.end(), &node1) == node0.neighbors.end())
	{
		return nullptr;
	}

	for (list<Edge2*>::const_iterator edges_iter = node0.inclInEdges.begin(), end = node0.inclInEdges.end();
		edges_iter != end; 
		edges_iter++)
	{
		if ((*edges_iter)->nodes[0] == &node0 && (*edges_iter)->nodes[1] == &node1 ||
			(*edges_iter)->nodes[0] == &node1 && (*edges_iter)->nodes[1] == &node0)
		{
			return *edges_iter;
		}
	}

	throw std::exception("Error in function Polycrystalline2::FindFreeEdge");
	return nullptr;
}

void Polycrystalline2::GenerateFreeSimplexesFromFreeNodes(size_t minNodesNumAxis_0)
{
	Simplex2* simp_buf = nullptr;
	
	for (size_t i = 0, max = _freeNodes.size() - minNodesNumAxis_0 - 1; i < max;)
	{
		for (size_t j = 0; j < minNodesNumAxis_0; j++)
		{
			_freeSimplexes.push_back(simp_buf = new Simplex2());
			simp_buf->edges[0] = new Edge2(*_freeNodes[i], *_freeNodes[i + 1]);
			simp_buf->edges[1] = new Edge2(*_freeNodes[i + 1], *_freeNodes[i + 1 + minNodesNumAxis_0]);
			simp_buf->edges[2] = new Edge2(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i]);
			_freeEdges.insert(
				_freeEdges.end(),
				{ simp_buf->edges[0],
				  simp_buf->edges[1],
				  simp_buf->edges[2] });

			simp_buf->edges[0]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[1]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[2]->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;

		if (!(i < max))
		{
			break;
		}

		for (size_t j = 0; j < minNodesNumAxis_0 - 1; j++)
		{
			_freeSimplexes.push_back(simp_buf = new Simplex2());
			simp_buf->edges[0] = new Edge2(*_freeNodes[i], *_freeNodes[i + 1]);
			simp_buf->edges[1] = new Edge2(*_freeNodes[i + 1], *_freeNodes[i + 1 + minNodesNumAxis_0]);
			simp_buf->edges[2] = new Edge2(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i]);
			_freeEdges.insert(
				_freeEdges.end(), 
				{ simp_buf->edges[0], 
				  simp_buf->edges[1], 
				  simp_buf->edges[2] });

			simp_buf->edges[0]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[1]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[2]->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;
	}

	Edge2* edge_buf = nullptr;
	
	for (size_t i = 1, max = _freeNodes.size() - minNodesNumAxis_0 - 1; i < max;)
	{
		for (size_t j = 0; j < minNodesNumAxis_0 - 1; j++)
		{
			_freeSimplexes.push_back(simp_buf = new Simplex2());
			if (edge_buf = FindFreeEdge(*_freeNodes[i], *_freeNodes[i + 1 + minNodesNumAxis_0]))
			{
				simp_buf->edges[0] = edge_buf;
			}
			else
			{
				simp_buf->edges[0] = new Edge2(*_freeNodes[i], *_freeNodes[i + 1 + minNodesNumAxis_0]);
				_freeEdges.push_back(simp_buf->edges[0]);
			}
			if (edge_buf = FindFreeEdge(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i + minNodesNumAxis_0]))
			{
				simp_buf->edges[1] = edge_buf;
			}
			else
			{
				simp_buf->edges[1] = new Edge2(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i + minNodesNumAxis_0]);
				_freeEdges.push_back(simp_buf->edges[1]);
			}
			if (edge_buf = FindFreeEdge(*_freeNodes[i], *_freeNodes[i + minNodesNumAxis_0]))
			{
				simp_buf->edges[2] = edge_buf;
			}
			else
			{
				simp_buf->edges[2] = new Edge2(*_freeNodes[i], *_freeNodes[i + minNodesNumAxis_0]);
				_freeEdges.push_back(simp_buf->edges[2]);
			}

			simp_buf->edges[0]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[1]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[2]->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;

		if (!(i < max))
		{
			break;
		}

		for (size_t j = 0; j < minNodesNumAxis_0; j++)
		{
			_freeSimplexes.push_back(simp_buf = new Simplex2());
			if (edge_buf = FindFreeEdge(*_freeNodes[i], *_freeNodes[i + 1 + minNodesNumAxis_0]))
			{
				simp_buf->edges[0] = edge_buf;
			}
			else
			{
				simp_buf->edges[0] = new Edge2(*_freeNodes[i], *_freeNodes[i + 1 + minNodesNumAxis_0]);
				_freeEdges.push_back(simp_buf->edges[0]);
			}
			if (edge_buf = FindFreeEdge(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i + minNodesNumAxis_0]))
			{
				simp_buf->edges[1] = edge_buf;
			}
			else
			{
				simp_buf->edges[1] = new Edge2(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i + minNodesNumAxis_0]);
				_freeEdges.push_back(simp_buf->edges[1]);
			}
			if (edge_buf = FindFreeEdge(*_freeNodes[i], *_freeNodes[i + minNodesNumAxis_0]))
			{
				simp_buf->edges[2] = edge_buf;
			}
			else
			{
				simp_buf->edges[2] = new Edge2(*_freeNodes[i], *_freeNodes[i + minNodesNumAxis_0]);
				_freeEdges.push_back(simp_buf->edges[2]);
			}

			simp_buf->edges[0]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[1]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[2]->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;
	}
}

void Polycrystalline2::GenerateFreeUniformMesh()
{
	double polycrysSizeAxis[2];
	polycrysSizeAxis[0] = _maxShellNodesCoor[0] - _minShellNodesCoor[0];
	polycrysSizeAxis[1] = _maxShellNodesCoor[1] - _minShellNodesCoor[1];

	size_t minNodesNumAxis[2];
	minNodesNumAxis[0] = (size_t)(polycrysSizeAxis[0] / _l_av) + 2;
	minNodesNumAxis[1] = (size_t)(polycrysSizeAxis[1] / (0.5 * sqrt(3) * _l_av)) + 2;

	GenerateFreeNodesEvenly(polycrysSizeAxis, minNodesNumAxis);
	GenerateFreeSimplexesFromFreeNodes(minNodesNumAxis[0]);
}

void Polycrystalline2::FitFreeNodesToShellNodes()
{
	size_t max = _shellNodes.size();
	//#pragma omp parallel for firstprivate(max)
	for (size_t i = 0; i < max; i++)
	{
		Node2* node_with_min_sqr_dist = _freeNodes[0];

		double min_sqr_dist = (*_freeNodes[0] - *_shellNodes[i]).SqrMagnitude();
		double sqr_magn_buf;
		for (size_t j = 1, max = _freeNodes.size(); j < max; j++)
		{
			if ((sqr_magn_buf = (*_freeNodes[j] - *_shellNodes[i]).SqrMagnitude()) < min_sqr_dist)
			{
				min_sqr_dist = sqr_magn_buf;
				node_with_min_sqr_dist = _freeNodes[j];
			}
		}

		//#pragma omp critical(FitFreeNodesToShellNodes)
		//{
		node_with_min_sqr_dist->SetPosition(_shellNodes[i]->GetPosition());
		node_with_min_sqr_dist->belongsToShellNode = _shellNodes[i];
		//}
	}
}

void Polycrystalline2::FitFreeNodesToShellEdges()
{
	size_t maxi = _freeEdges.size();
	//#pragma omp parallel for firstprivate(maxi)
	for (size_t i = 0; i < maxi; i++)
	{
		for (size_t j = 0, maxj = _shellEdges.size(); j < maxj; j++)
		{
			if (_freeEdges[i]->nodes[0]->belongsToShellNode ||
				_freeEdges[i]->nodes[1]->belongsToShellNode)
			{
				continue;
			}
			if (((*_freeEdges[i]->nodes[0])[0] < (*_shellEdges[j]->nodes[0])[0]  &&
				 (*_freeEdges[i]->nodes[1])[0] < (*_shellEdges[j]->nodes[0])[0]  &&
				 (*_freeEdges[i]->nodes[0])[0] < (*_shellEdges[j]->nodes[1])[0]  &&
				 (*_freeEdges[i]->nodes[1])[0] < (*_shellEdges[j]->nodes[1])[0]) ||

				((*_freeEdges[i]->nodes[0])[1] < (*_shellEdges[j]->nodes[0])[1]  &&
				 (*_freeEdges[i]->nodes[1])[1] < (*_shellEdges[j]->nodes[0])[1]  &&
				 (*_freeEdges[i]->nodes[0])[1] < (*_shellEdges[j]->nodes[1])[1]  &&
				 (*_freeEdges[i]->nodes[1])[1] < (*_shellEdges[j]->nodes[1])[1]) ||

				((*_freeEdges[i]->nodes[0])[0] > (*_shellEdges[j]->nodes[0])[0]  &&
				 (*_freeEdges[i]->nodes[1])[0] > (*_shellEdges[j]->nodes[0])[0]  &&
				 (*_freeEdges[i]->nodes[0])[0] > (*_shellEdges[j]->nodes[1])[0]  &&
				 (*_freeEdges[i]->nodes[1])[0] > (*_shellEdges[j]->nodes[1])[0]) ||

				((*_freeEdges[i]->nodes[0])[1] > (*_shellEdges[j]->nodes[0])[1]  &&
				 (*_freeEdges[i]->nodes[1])[1] > (*_shellEdges[j]->nodes[0])[1]  &&
				 (*_freeEdges[i]->nodes[0])[1] > (*_shellEdges[j]->nodes[1])[1]  &&
				 (*_freeEdges[i]->nodes[1])[1] > (*_shellEdges[j]->nodes[1])[1]))
			{
				continue;
			}
			//_freeEdges[i]->nodes[0]->SetPosition(0, 0);
			//_freeEdges[i]->nodes[1]->SetPosition(0, 0);
			double A1 = (*_freeEdges [i]->nodes[0])[1] - (*_freeEdges [i]->nodes[1])[1];
			double B1 = (*_freeEdges [i]->nodes[1])[0] - (*_freeEdges [i]->nodes[0])[0]; 
			double C1 = -A1 * (*_freeEdges [i]->nodes[0])[0] - B1 * (*_freeEdges [i]->nodes[0])[1];
			double A2 = (*_shellEdges[j]->nodes[0])[1] - (*_shellEdges[j]->nodes[1])[1];
			double B2 = (*_shellEdges[j]->nodes[1])[0] - (*_shellEdges[j]->nodes[0])[0];
			double C2 = -A2 * (*_shellEdges[j]->nodes[0])[0] - B2 * (*_shellEdges[j]->nodes[0])[1];
			double minus_inv_zn = -1.0 / DET(A1, B1, A2, B2);
			double x = DET(C1, B1, C2, B2) * minus_inv_zn;
			double y = DET(A1, C1, A2, C2) * minus_inv_zn;
			if ((std::min((*_freeEdges[i]->nodes[0])[0], (*_freeEdges[i]->nodes[1])[0]) - DBL_EPSILON < x && x < std::max((*_freeEdges[i]->nodes[0])[0], (*_freeEdges[i]->nodes[1])[0]) + DBL_EPSILON) &&//BETWEEN((*_freeEdges [i]->nodes[0])[0], (*_freeEdges [i]->nodes[1])[0], x) &&
				BETWEEN((*_freeEdges [i]->nodes[0])[1], (*_freeEdges [i]->nodes[1])[1], y) &&
				BETWEEN((*_shellEdges[j]->nodes[0])[0], (*_shellEdges[j]->nodes[1])[0], x) &&
				BETWEEN((*_shellEdges[j]->nodes[0])[1], (*_shellEdges[j]->nodes[1])[1], y))
			{
				Vector2 intersect_point(x, y);
				if ((_freeEdges[i]->nodes[0]->GetPosition() - intersect_point).SqrMagnitude() <
					(_freeEdges[i]->nodes[1]->GetPosition() - intersect_point).SqrMagnitude())
				{
					//#pragma omp critical(FitFreeNodesToShellEdges)
					//{
					_freeEdges[i]->nodes[0]->SetPosition(intersect_point);
					_freeEdges[i]->nodes[0]->belongsToShellEdge = _shellEdges[j];
					//}
				}
				else
				{
					//#pragma omp critical(FitFreeNodesToShellEdges)
					//{
					_freeEdges[i]->nodes[1]->SetPosition(intersect_point);
					_freeEdges[i]->nodes[1]->belongsToShellEdge = _shellEdges[j];
					//}
				}
			}
		}
	}
}

void Polycrystalline2::FitFreeMeshToShells()
{
	FitFreeNodesToShellNodes();
	FitFreeNodesToShellEdges();
}

const bool Polycrystalline2::Contains(const Node2& node) const
{
	for (auto crys : crystallites)
	{
		if (crys->Contains(node))
		{
			return true;
		}
	}

	return false;
}

void AddNodesData(ofstream& nodesData, const list<Crystallite2*>& crystallites)
{
	// Can be parallelized, but it's a bad idea.
	size_t index = 0;
	Vector2 bufPos;
	for (list<Crystallite2*>::const_iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		for (list<Simplex2*>::const_iterator simp_iter = (*crys_iter)->simplexes.begin();
			simp_iter != (*crys_iter)->simplexes.end();
			simp_iter++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (!(*simp_iter)->edges[i]->nodes[j]->isAddedToNodesData)
					{
						bufPos = (*simp_iter)->edges[i]->nodes[j]->GetPosition();
						if (crys_iter != crystallites.begin() ||
							simp_iter != (*crys_iter)->simplexes.begin() ||
							i > 0 ||
							j > 0)
						{
							nodesData << '\n';
						}
						nodesData << bufPos[0] << ' ';
						nodesData << bufPos[1];
						(*simp_iter)->edges[i]->nodes[j]->isAddedToNodesData = true;
						(*simp_iter)->edges[i]->nodes[j]->globalNum = index;
						index++;
					}
				}
			}
		}
	}
}

void AddNodesData(vector<double>& nodesData, const list<Crystallite2*>& crystallites)
{
	// Can be parallelized, but it's a bad idea.
	size_t index = 0;
	Vector2 bufPos;
	for (list<Crystallite2*>::const_iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		for (list<Simplex2*>::const_iterator simp_iter = (*crys_iter)->simplexes.begin();
			simp_iter != (*crys_iter)->simplexes.end();
			simp_iter++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (!(*simp_iter)->edges[i]->nodes[j]->isAddedToNodesData)
					{
						bufPos = (*simp_iter)->edges[i]->nodes[j]->GetPosition();
						nodesData.push_back(bufPos[0]);
						nodesData.push_back(bufPos[1]);
						(*simp_iter)->edges[i]->nodes[j]->isAddedToNodesData = true;
						(*simp_iter)->edges[i]->nodes[j]->globalNum = index;
						index++;
					}
				}
			}
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
	if (simp.edges[1]->nodes[0] == simp.edges[0]->nodes[0] ||
		simp.edges[1]->nodes[0] == simp.edges[0]->nodes[1])
	{
		fenData << simp.edges[1]->nodes[1]->globalNum;
	}
	else
	{
		fenData << simp.edges[1]->nodes[0]->globalNum;
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

void AddFENodesData(ofstream& feNodesData, const list<Crystallite2*>& crystallites)
{
	// Can be parallelized, but it's a bad idea.
	bool isFirst = true;
	for (list<Crystallite2*>::const_iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		size_t simp_size = (*crys_iter)->simplexes.size();
		if (*crys_iter && (*crys_iter)->simplexes.size() > 0)
		{
			for (list<unique_ptr<Simplex2>*>::const_iterator simp_iter = (*crys_iter)->simplexes.begin();
				simp_iter != (*crys_iter)->simplexes.end();
				simp_iter++)
			{
				if (*simp_iter)
				{
					AddFENodesDataFromSimpex(feNodesData, ***simp_iter, isFirst);
					isFirst = false;
				}
			}
		}
	}
}

void AddFENodesData(vector<size_t>& feNodesData, const list<Crystallite2*>& crystallites)
{
	// Can be parallelized, but it's a bad idea.
	for (list<Crystallite2*>::const_iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		if (*crys_iter)
		{
			for (list<unique_ptr<Simplex2>*>::const_iterator simp_iter = (*crys_iter)->simplexes.begin();
				simp_iter != (*crys_iter)->simplexes.end();
				simp_iter++)
			{
				if (*simp_iter)
				{
					AddFENodesDataFromSimpex(feNodesData, ***simp_iter);
				}
			}
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
			new ShellNode2(
				dbuf[0], 
				dbuf[1])->GetPtrToUniquePtr());
	}
	shell_nodes.clear();
	shell_nodes.seekg(0);

	size_t ibuf[2];
	size_t numsNum;
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
		cryses_edges >> numsNum;
		//numsNum *= 2; // numsNum *= 3; for 3 dimensions.
		for (size_t j = 0; j < numsNum; j++)
		{
			cryses_edges >> ibuf[0];
			cryses_edges >> ibuf[1];
			
			for (auto &s_edge : (*_shellNodes[ibuf[0]])->inclInEdges)
			{
				//if (s_edge->nodes[1] == _shellNodes[ibuf[1]] ||
				//	  s_edge->nodes[1] == _shellNodes[ibuf[0]])
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
					new unique_ptr<ShellEdge2>(
						new ShellEdge2(
							**_shellNodes[ibuf[0]],
							**_shellNodes[ibuf[1]])));
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
			new ShellNode2(
				shell_nodes[i],
				shell_nodes[i + 1])->GetPtrToUniquePtr());

		size_t ibuf[2];
		size_t numsNum;
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
			numsNum = *(cp_iter++); // numsNum = *cp_iter; cp_iter++;
			//numsNum *= 2; // numsNum *= 3; for 3 dimensions.
			for (size_t j = 0; j < numsNum; j++)
			{
				for (auto &s_edge : (*_shellNodes[ibuf[0]])->inclInEdges)
				{
					//if (s_edge->nodes[1] == _shellNodes[ibuf[1]] ||
					//	  s_edge->nodes[0] == _shellNodes[ibuf[1]])
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
							new ShellEdge2(
								**_shellNodes[*(cp_iter)],
								**_shellNodes[*(cp_iter + 1)])->GetPtrToUniquePtr());
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
	AddNodesData(nodesData, crystallites);
	AddFENodesData(feNodesData, crystallites);
}

void Polycrystalline2::OutputData(vector<double>& nodesData, vector<size_t>& feNodesData) const
{
	AddNodesData(nodesData, crystallites);
	AddFENodesData(feNodesData, crystallites);
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
		}
	}
}