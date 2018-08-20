#include <iostream>
#include <fstream>
#include "Polycrystalline2.h"
using namespace std;

int main(int argc, char* argv[])
{
	std::ifstream cryses_nodes_pos("Input/cryses_nodes_pos.txt");
	std::ifstream cryses_edges("Input/cryses_edges.txt");
	std::ifstream gener_params("Input/gener_params.txt");

	Polycrystalline2 polycr(cryses_nodes_pos, cryses_edges, gener_params);

	cryses_nodes_pos.close();
	cryses_edges.close();
	gener_params.close();

	polycr.GenerateFreeUniformMesh();
	polycr.FitFreeMeshToShells();
	//polycr.DeleteExternalNodes();
	int input;
	while (true)
	{
		cout << "\n1. Distribute nodes evenly.\n"
				"2. Make simplexes equilateral.\n"
				"3. Divide extended simplexes.\n"
				"4. Close.\n";
		cin >> input;
		if (input == 1)
		{
			polycr.DistributeNodesEvenly();
		}
		else if (input == 2)
		{
			polycr.MakeSimplexesEquilateral();
		}
		else if (input == 3)
		{
			polycr.DivideExtendedSimplexes();
		}
		else
		{
			break;
		}
	}

	std::ofstream nodes_pos("Output/nodes_pos.txt");
	std::ofstream fe_nodes("Output/fe_nodes.txt");

	polycr.OutputData(nodes_pos, fe_nodes);

	nodes_pos.close();
	fe_nodes.close();


	return 0;
}