#include <iostream>
#include <fstream>
#include "Polycrystalline3.h"
using namespace std;

int main(int argc, char* argv[])
{
	Polycrystalline3 polycr;
	polycr.InputData();
	polycr.TriangulateShell();
	polycr.SetLinksWithShell();
	polycr.DelaunayPostprocessing();
	polycr.TriangulateCrystallitesVolumes();
	polycr.OutputData();

	return 0;
}