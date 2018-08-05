#include <iostream>
#include <fstream>
#include "Polycrystalline2.h"

using namespace std;


int main(int argc, char* argv[])
{
	ifstream cryses_nodes_pos("Input/cryses_nodes_pos.txt");
	ifstream cryses_edges("Input/cryses_edges.txt");
	ifstream gener_params("Input/gener_params.txt");



	cryses_nodes_pos.close();
	cryses_edges.close();
	gener_params.close();

	ofstream nodes_pos("Output/nodes_pos.txt");
	ofstream fe_nodes("Output/fe_nodes.txt");



	nodes_pos.close();
	fe_nodes.close();

	cout << '\n';
	system("pause");
	return 0;
}