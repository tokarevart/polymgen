#include <iostream>
#include <fstream>
#include "Polycrystalline2.h"

using namespace std;


int main(int argc, char* argv[])
{
	ifstream in_nodes("Input/nodes.txt");
	ifstream crystallites("Input/crystallites.txt");
	ifstream crys_edges("Input/crys_edges.txt");


	in_nodes.close();
	crystallites.close();
	crys_edges.close();

	cout << '\n';
	system("pause");
	return 0;
}