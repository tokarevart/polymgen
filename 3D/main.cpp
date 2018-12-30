#include <iostream>
#include <fstream>
#include "Polycrystal3.h"
using namespace std;

int main(int argc, char* argv[])
{
	Polycrystal3 polycr("input_cube_test.txt");
	polycr.generateMeshNoStruct(0.16);
	polycr.outputData();
	return 0;
}