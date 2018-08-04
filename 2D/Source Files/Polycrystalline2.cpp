#include "Polycrystalline2.h"

using namespace std;


void Polycrystalline2::AddCrystallite()
{
	Crystallite2* buf;
	crystallites.push_back(buf = new Crystallite2());
	buf->inclInPolycrys = this;
}

void Polycrystalline2::AddCrystallite(Crystallite2* const& crys)
{
	crystallites.push_back(crys);
	crys->inclInPolycrys = this;
}

const bool Polycrystalline2::IsContaining(const Crystallite2& crys)
{
	return find(crystallites.begin(), crystallites.end(), &crys) != crystallites.end();
}

const bool Polycrystalline2::IsContaining(const Simplex2& simp)
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

const bool Polycrystalline2::IsContaining(const Edge2& edge)
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

const bool Polycrystalline2::IsContaining(const Node2& node)
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

const bool Polycrystalline2::Contains(const Simplex2& simp)
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

const bool Polycrystalline2::Contains(const Edge2& edge)
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

const bool Polycrystalline2::Contains(const Node2& node)
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
void Polycrystalline2::InputData(ifstream& shell_nodes, ifstream& crys_edges)
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
	for (list<Crystallite2*>::iterator crys_iter; !crys_edges.eof(); crys_iter++)
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
		crys_edges >> numsNum;
		numsNum *= 2; // numsNum *= 3; for 3 dimensions.
		for (size_t j = 0; j < numsNum; j++)
		{
			crys_edges >> ibuf[0];
			crys_edges >> ibuf[1];
			(*crys_iter)->shellEdges.push_back(
				new ShellEdge2(
					*v_shell_nodes[ibuf[0]], 
					*v_shell_nodes[ibuf[1]]));
		}
	}
	crys_edges.clear();
	crys_edges.seekg(0);
}

// Different for 3 dimensions.
void Polycrystalline2::InputData(const vector<double>& shell_nodes, const vector<size_t>& crys_edges)
{
	vector<ShellNode2*> v_shell_nodes;

	for (size_t i = 0, max = shell_nodes.size(); i < max; i += 2) // i += 3 for 3 dimensions
	{
		v_shell_nodes.push_back(
			new ShellNode2(
				shell_nodes[i], 
				shell_nodes[i + 1]));

		size_t ibuf[2];
		size_t numsNum;
		vector<size_t>::const_iterator cp_iter = crys_edges.begin();
		for (list<Crystallite2*>::iterator crys_iter; cp_iter != crys_edges.end(); crys_iter++)
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
}

void Polycrystalline2::OutputData(ofstream& nodesData, ofstream& feNodesData)
{
	AddNodesData(nodesData, crystallites);
	AddFENodesData(feNodesData, crystallites);
}

void Polycrystalline2::OutputData(vector<double>& nodesData, vector<size_t>& feNodesData)
{
	AddNodesData(nodesData, crystallites);
	AddFENodesData(feNodesData, crystallites);
}

Polycrystalline2::Polycrystalline2()
{
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