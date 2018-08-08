#include <iostream>
#include <fstream>
#include "Polycrystalline2.h"
#include "Timer.h"


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
	polycr.FitFreeNodesToShellNodes();
	polycr.Debug();

	std::ofstream nodes_pos("Output/nodes_pos.txt");
	std::ofstream fe_nodes("Output/fe_nodes.txt");

	polycr.OutputData(nodes_pos, fe_nodes);

	nodes_pos.close();
	fe_nodes.close();

	//std::cout << '\n';
	//system("pause");
	return 0;
}