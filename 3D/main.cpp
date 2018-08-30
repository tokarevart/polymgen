#include <iostream>
#include <fstream>
#include "Polycrystalline3.h"
using namespace std;

int main(int argc, char* argv[])
{
	Polycrystalline3 polycr;
	polycr.InputData();
	polycr.TriangulateShell();
	polycr.DistributeVertexesOnShellEvenly(10);
	polycr.DelaunayPostprocessing();
	polycr.DistributeVertexesOnShellEvenly(10);
	polycr.OutputData();

	return 0;
}