#include <iostream>
#include <fstream>
#include "Polycrystal3.h"

int main(int argc, char* argv[])
{
	double preferredEdgeLength;
	std::cout << "Enter preferred edge length: ";
	std::cin >> preferredEdgeLength;

	std::cout << "Generating mesh...";
	Polycrystal3 polycr("input_cube_test.txt");
	polycr.generateMeshNoStruct(preferredEdgeLength); // 0.45 is quite good
	polycr.outputData();
	return 0;
}