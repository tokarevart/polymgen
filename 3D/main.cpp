#include <iostream>
#include <fstream>
#include "Polycrystal3.h"
using namespace std;

int main(int argc, char* argv[])
{
	Polycrystal3 polycr("input_cube_test.txt");
	polycr.TriangulatePolycrystalNoStruct(0.45);
	polycr.OutputData();
	return 0;
}