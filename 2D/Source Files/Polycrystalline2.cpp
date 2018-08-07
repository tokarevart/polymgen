#include "Polycrystalline2.h"

using namespace std;


size_t Polycrystalline2::GenerateFreeNodesEvenly(double* const polycrysSizeAxis, size_t* const minNodesNumAxis)
{
	double startGenCoor[2];
	startGenCoor[0] = _minShellNodesCoor[0] + (polycrysSizeAxis[0] - (minNodesNumAxis[0] + 1) * _l_av) * 0.5;
	startGenCoor[1] = _minShellNodesCoor[1] + (polycrysSizeAxis[1] - minNodesNumAxis[1] * _l_av) * 0.5;

	size_t freeNodesNum = minNodesNumAxis[0] * minNodesNumAxis[1] + (minNodesNumAxis[1] + 1) / 2;

	double node_coors[2];
	node_coors[0] = startGenCoor[0];
	node_coors[1] = startGenCoor[1];
	_freeNodes.reserve(freeNodesNum);
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

		if (!i < freeNodesNum)
		{
			break;
		}

		node_coors[0] -= minNodesNumAxis[0] * _l_av + _l_av * 0.5;
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
	}

	return freeNodesNum;
}

Edge2* Polycrystalline2::FindFreeEdge(const Node2& node0, const Node2& node1)
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

	throw exception("Error in function Polycrystalline2::FindFreeEdge");
	return nullptr;
}

void Polycrystalline2::GenerateFreeSimplexesFromFreeNodes(size_t minNodesNumAxis_0)
{
	Simplex2* simp_buf;
	
	for (size_t i = 0, max = _freeNodes.size() - minNodesNumAxis_0 - 1; i < max;)
	{
		for (size_t j = 0; j < minNodesNumAxis_0; j++)
		{
			_freeSimplexes.push_back(simp_buf = new Simplex2());
			simp_buf->edges[0] = new Edge2(*_freeNodes[i], *_freeNodes[i + 1]);
			simp_buf->edges[1] = new Edge2(*_freeNodes[i + 1], *_freeNodes[i + 1 + minNodesNumAxis_0]);
			simp_buf->edges[2] = new Edge2(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i]);

			simp_buf->edges[0]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[1]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[2]->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;

		if (!i < max)
		{
			break;
		}

		for (size_t j = 0; j < minNodesNumAxis_0 - 1; j++)
		{
			_freeSimplexes.push_back(simp_buf = new Simplex2());
			simp_buf->edges[0] = new Edge2(*_freeNodes[i], *_freeNodes[i + 1]);
			simp_buf->edges[1] = new Edge2(*_freeNodes[i + 1], *_freeNodes[i + 1 + minNodesNumAxis_0]);
			simp_buf->edges[2] = new Edge2(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i]);

			simp_buf->edges[0]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[1]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[2]->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;
	}

	Edge2* edge_buf;
	
	for (size_t i = 0, max = _freeNodes.size() - minNodesNumAxis_0 - 1; i < max;)
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
			}
			if (edge_buf = FindFreeEdge(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i + minNodesNumAxis_0]))
			{
				simp_buf->edges[1] = edge_buf;
			}
			else
			{
				simp_buf->edges[1] = new Edge2(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i + minNodesNumAxis_0]);
			}
			if (edge_buf = FindFreeEdge(*_freeNodes[i], *_freeNodes[i + minNodesNumAxis_0]))
			{
				simp_buf->edges[2] = edge_buf;
			}
			else
			{
				simp_buf->edges[2] = new Edge2(*_freeNodes[i], *_freeNodes[i + minNodesNumAxis_0]);
			}

			simp_buf->edges[0]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[1]->inclInSimplexes.push_back(simp_buf);
			simp_buf->edges[2]->inclInSimplexes.push_back(simp_buf);

			i++;
		}
		i++;

		if (!i < max)
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
			}
			if (edge_buf = FindFreeEdge(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i + minNodesNumAxis_0]))
			{
				simp_buf->edges[1] = edge_buf;
			}
			else
			{
				simp_buf->edges[1] = new Edge2(*_freeNodes[i + 1 + minNodesNumAxis_0], *_freeNodes[i + minNodesNumAxis_0]);
			}
			if (edge_buf = FindFreeEdge(*_freeNodes[i], *_freeNodes[i + minNodesNumAxis_0]))
			{
				simp_buf->edges[2] = edge_buf;
			}
			else
			{
				simp_buf->edges[2] = new Edge2(*_freeNodes[i], *_freeNodes[i + minNodesNumAxis_0]);
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
	minNodesNumAxis[0] = (size_t)(polycrysSizeAxis[0] / _l_av) + 1;
	minNodesNumAxis[1] = (size_t)(polycrysSizeAxis[1] / _l_av) + 1;

	size_t freeNodesNum = GenerateFreeNodesEvenly(polycrysSizeAxis, minNodesNumAxis);

	_freeEdges.reserve(3 + 2 * (_freeNodes.size() - 3) + minNodesNumAxis[1] - 2);
	GenerateFreeSimplexesFromFreeNodes(minNodesNumAxis[0]);
}

const bool Polycrystalline2::IsContaining(const Crystallite2& crys) const
{
	return find(crystallites.begin(), crystallites.end(), &crys) != crystallites.end();
}

const bool Polycrystalline2::IsContaining(const Simplex2& simp) const
{
	for (auto crys : crystallites)
	{
		if (crys->IsContaining(simp))
		{
			return true;
		}
	}

	return false;
}

const bool Polycrystalline2::IsContaining(const Edge2& edge) const
{
	for (auto crys : crystallites)
	{
		if (crys->IsContaining(edge))
		{
			return true;
		}
	}

	return false;
}

const bool Polycrystalline2::IsContaining(const Node2& node) const
{
	for (auto crys : crystallites)
	{
		if (crys->IsContaining(node))
		{
			return true;
		}
	}

	return false;
}

const bool Polycrystalline2::Contains(const Simplex2& simp) const
{
	for (auto crys : crystallites)
	{
		if (crys->Contains(simp))
		{
			return true;
		}
	}

	return false;
}

const bool Polycrystalline2::Contains(const Edge2& edge) const
{
	for (auto crys : crystallites)
	{
		if (crys->Contains(edge))
		{
			return true;
		}
	}

	return false;
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
						nodesData << bufPos[0] << ' ';
						nodesData << bufPos[1] << '\n';
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
void AddFENodesDataFromSimpex(ofstream& fenData, const Simplex2& simp)
{
	fenData << simp.edges[0]->nodes[0]->globalNum << ' ';
	fenData << simp.edges[0]->nodes[1]->globalNum << ' ';
	if (simp.edges[1]->nodes[0] == simp.edges[0]->nodes[0] ||
		simp.edges[1]->nodes[0] == simp.edges[0]->nodes[1])
	{
		fenData << simp.edges[1]->nodes[1]->globalNum;
	}
	else
	{
		fenData << simp.edges[1]->nodes[0]->globalNum;
	}
	fenData << '\n';
}

// Different for 3 dimensions.
void AddFENodesDataFromSimpex(vector<size_t>& fenData, const Simplex2& simp)
{
	fenData.push_back(simp.edges[0]->nodes[0]->globalNum);
	fenData.push_back(simp.edges[0]->nodes[1]->globalNum);
	if (simp.edges[1]->nodes[0] == simp.edges[0]->nodes[0] ||
		simp.edges[1]->nodes[0] == simp.edges[0]->nodes[1])
	{
		fenData.push_back(simp.edges[1]->nodes[1]->globalNum);
	}
	else
	{
		fenData.push_back(simp.edges[1]->nodes[0]->globalNum);
	}
}

void AddFENodesData(ofstream& feNodesData, const list<Crystallite2*>& crystallites)
{
	// Can be parallelized, but it's a bad idea.
	for (list<Crystallite2*>::const_iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		if (*crys_iter)
		{
			for (list<Simplex2*>::const_iterator simp_iter = (*crys_iter)->simplexes.begin();
				simp_iter != (*crys_iter)->simplexes.end();
				simp_iter++)
			{
				if (*simp_iter)
				{
					AddFENodesDataFromSimpex(feNodesData, **simp_iter);
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
			for (list<Simplex2*>::const_iterator simp_iter = (*crys_iter)->simplexes.begin();
				simp_iter != (*crys_iter)->simplexes.end();
				simp_iter++)
			{
				if (*simp_iter)
				{
					AddFENodesDataFromSimpex(feNodesData, **simp_iter);
				}
			}
		}
	}
}

// Different for 3 dimensions.
void Polycrystalline2::InputData(ifstream& shell_nodes, ifstream& cryses_edges, ifstream& gener_params)
{
	vector<ShellNode2*> v_shell_nodes;

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
		v_shell_nodes.push_back(
			new ShellNode2(
				dbuf[0], 
				dbuf[1]));
	}
	shell_nodes.clear();
	shell_nodes.seekg(0);

	size_t ibuf[2];
	size_t numsNum;
	for (list<Crystallite2*>::iterator crys_iter; !cryses_edges.eof(); crys_iter++)
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
		numsNum *= 2; // numsNum *= 3; for 3 dimensions.
		for (size_t j = 0; j < numsNum; j++)
		{
			cryses_edges >> ibuf[0];
			cryses_edges >> ibuf[1];
			(*crys_iter)->shellEdges.push_back(
				new ShellEdge2(
					*v_shell_nodes[ibuf[0]], 
					*v_shell_nodes[ibuf[1]]));
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
	vector<ShellNode2*> v_shell_nodes;

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

		v_shell_nodes.push_back(
			new ShellNode2(
				shell_nodes[i], 
				shell_nodes[i + 1]));

		size_t ibuf[2];
		size_t numsNum;
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
			numsNum *= 2; // numsNum *= 3; for 3 dimensions.
			for (size_t j = 0; j < numsNum; j++)
			{
				(*crys_iter)->shellEdges.push_back(
					new ShellEdge2(
						*v_shell_nodes[*(cp_iter)],
						*v_shell_nodes[*(cp_iter + 1)]));

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