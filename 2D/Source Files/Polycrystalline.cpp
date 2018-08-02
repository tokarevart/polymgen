#include "Polycrystalline.h"

using namespace std;


void Polycrystalline::AddCrystallite()
{
	Crystallite2* buf;
	crystallites.push_back(buf = new Crystallite2());
	buf->inclInPolycrys = this;
}

void Polycrystalline::AddCrystallite(Crystallite2* &crys)
{
	crystallites.push_back(crys);
	crys->inclInPolycrys = this;
}

bool Polycrystalline::IsContaining(Simplex2 &simp)
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

bool Polycrystalline::IsContaining(Edge &edge)
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

bool Polycrystalline::IsContaining(Node2 &node)
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

bool Polycrystalline::Contains(Crystallite2 &crys)
{
	return find(crystallites.begin(), crystallites.end(), &crys) != crystallites.end();
}

bool Polycrystalline::Contains(Simplex2 &simp)
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

bool Polycrystalline::Contains(Edge &edge)
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

bool Polycrystalline::Contains(Node2 &node)
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

void AddNodesData(ofstream &nodesData, list<Crystallite2*> &crystallites)
{
	// Can be parallelized, but it's a bad idea.
	int index = 0;
	for (list<Crystallite2*>::iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		for (list<Simplex2*>::iterator simp_iter = (*crys_iter)->simplexes.begin();
			simp_iter != (*crys_iter)->simplexes.end();
			simp_iter++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (!(*simp_iter)->edges[i]->nodes[j]->isAddedToNodesData)
					{
						nodesData << (*simp_iter)->edges[i]->nodes[j]->coors[0] << ' ';
						nodesData << (*simp_iter)->edges[i]->nodes[j]->coors[1] << '\n';
						(*simp_iter)->edges[i]->nodes[j]->isAddedToNodesData = true;
						(*simp_iter)->edges[i]->nodes[j]->globalNum = index;
						index++;
					}
				}
			}
		}
	}
}

void AddNodesData(vector<double> &nodesData, list<Crystallite2*> &crystallites)
{
	// Can be parallelized, but it's a bad idea.
	int index = 0;
	for (list<Crystallite2*>::iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		for (list<Simplex2*>::iterator simp_iter = (*crys_iter)->simplexes.begin();
			simp_iter != (*crys_iter)->simplexes.end();
			simp_iter++)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 2; j++)
				{
					if (!(*simp_iter)->edges[i]->nodes[j]->isAddedToNodesData)
					{
						nodesData.push_back((*simp_iter)->edges[i]->nodes[j]->coors[0]);
						nodesData.push_back((*simp_iter)->edges[i]->nodes[j]->coors[1]);
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
void AddFENodesDataFromSimpex(ofstream &fenData, Simplex2 &simp)
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
void AddFENodesDataFromSimpex(vector<size_t> &fenData, Simplex2 &simp)
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

void AddFENodesData(ofstream &feNodesData, list<Crystallite2*> &crystallites)
{
	// Can be parallelized, but it's a bad idea.
	for (list<Crystallite2*>::iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		for (list<Simplex2*>::iterator simp_iter = (*crys_iter)->simplexes.begin();
			simp_iter != (*crys_iter)->simplexes.end();
			simp_iter++)
		{
			AddFENodesDataFromSimpex(feNodesData, **simp_iter);
		}
	}
}

void AddFENodesData(vector<size_t> &feNodesData, list<Crystallite2*> &crystallites)
{
	// Can be parallelized, but it's a bad idea.
	for (list<Crystallite2*>::iterator crys_iter = crystallites.begin();
		crys_iter != crystallites.end();
		crys_iter++)
	{
		for (list<Simplex2*>::iterator simp_iter = (*crys_iter)->simplexes.begin();
			simp_iter != (*crys_iter)->simplexes.end();
			simp_iter++)
		{
			AddFENodesDataFromSimpex(feNodesData, **simp_iter);
		}
	}
}

// Different for 3 dimensions.
void Polycrystalline::InputData(ifstream &shell_nodes, ifstream &crys_planes)
{
	vector<ShellNode2*> v_shell_nodes;

	double dbuf[2];
	while (!shell_nodes.eof())
	{
		shell_nodes >> dbuf[0];
		shell_nodes >> dbuf[1];
		v_shell_nodes.push_back(
			new ShellNode2(
				dbuf[0], 
				dbuf[1]));
	}
	shell_nodes.clear();
	shell_nodes.seekg(0);

	size_t ibuf[2];
	size_t numsNum;
	for (list<Crystallite2*>::iterator crys_iter; !crys_planes.eof(); crys_iter++)
	{
		if (crystallites.empty())
		{
			// crystallites.push_back(new Crystallite2()); instead of AddCrystallite();
			AddCrystallite();
			crys_iter = crystallites.begin();
		}
		else
		{
			// crystallites.push_back(new Crystallite2()); instead of AddCrystallite();
			AddCrystallite();
			crys_iter++;
		}
		crys_planes >> numsNum;
		numsNum *= 2; // numsNum *= 3; for 3 dimensions.
		for (int j = 0; j < numsNum; j++)
		{
			crys_planes >> ibuf[0];
			crys_planes >> ibuf[1];
			(*crys_iter)->shellEdges.push_back(
				new ShellEdge(
					*v_shell_nodes[ibuf[0]], 
					*v_shell_nodes[ibuf[1]]));
		}
	}
}

// Different for 3 dimensions.
void Polycrystalline::InputData(vector<double> &shell_nodes, vector<size_t> &crys_planes)
{
}

void Polycrystalline::OutputData(ofstream &nodesData, ofstream &feNodesData)
{
}

void Polycrystalline::OutputData(vector<double> &nodesData, vector<size_t> &feNodesData)
{
	AddNodesData(nodesData, crystallites);
	AddFENodesData(feNodesData, crystallites);
}

Polycrystalline::Polycrystalline()
{
}

Polycrystalline::~Polycrystalline()
{
	for (auto &crys : crystallites)
	{
		if (crys)
		{
			delete crys;
		}
	}
}