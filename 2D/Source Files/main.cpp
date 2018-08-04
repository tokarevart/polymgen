#include <iostream>
#include <fstream>
#include "Polycrystalline2.h"

using namespace std;


int main(int argc, char* argv[])
{
	ifstream in_nodes("Input/nodes.txt");
	ifstream crys_edges("Input/crys_edges.txt");
	ifstream gener_params("Input/gener_params.txt");



	in_nodes.close();
	crys_edges.close();
	gener_params.close();

	cout << '\n';
	system("pause");
	return 0;
}