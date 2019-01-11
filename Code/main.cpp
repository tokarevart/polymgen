#include <iostream>
#include <fstream>
#include "Polycrystal3.h"

int main(int argc, char* argv[])
{
	double preferredEdgeLength;
	std::cout << "Enter preferred edge length: ";
	std::cin >> preferredEdgeLength;

	std::cout << "\nInitializing polycrystal data...";
	Polycrystal3 polycr("input_cube_test.txt");
	std::cout << " done.\n";

	std::cout << "Generating mesh...";
	PolycrMesh* mesh = polycr.generateMesh(preferredEdgeLength); // 0.45 is quite good
	std::cout << " done.\n";

	std::cout << "Outputting data to file...";
	polycr.outputData();
	std::cout << " done.\n";

	std::ofstream polycr_mesh("polycr_mesh_test.txt");
	polycr_mesh << mesh->nodesNum << '\n';
	for (size_t i = 0ull; i < 3ull * mesh->nodesNum; i += 3ull)
		polycr_mesh
			<< mesh->nodesPositions[i] << ' '
			<< mesh->nodesPositions[i + 1ull] << ' '
			<< mesh->nodesPositions[i + 2ull] << '\n';
	polycr_mesh << '\n' << mesh->tetrsNum << '\n';
	polycr_mesh << mesh->crysesNum << '\n';
	for (size_t i = 0ull; i < mesh->crysesNum; i++)
		polycr_mesh << mesh->crysesTetrsNum[i] << ' ';
	polycr_mesh << "\n...";
	delete mesh;
	return 0;
}