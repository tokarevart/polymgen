#include <iostream>
#include <fstream>
#include "Polycrystal3.h"
#include "PolycrGen.h"

int main(int argc, char* argv[])
{
    double preferredEdgeLength = 0.45;
    //std::cout << "Enter preferred edge length: ";
    //std::cin >> preferredEdgeLength;

    #define GENERATE
    #ifdef GENERATE

    std::cout << "Generating polycrystal...";
    CrysesShell* shell = polycrgen::generateCuboidsPolycrystal(9, 9, 9);
    std::cout << " done.\n";

    //std::ofstream polycr_mesh("polycr_mesh_test.txt");
    //polycr_mesh << shell->nNodes << '\n';
    //for (size_t i = 0ull; i < 3ull * shell->nNodes; i += 3ull)
    //    polycr_mesh
    //        << shell->nodesPositions[i] << ' '
    //        << shell->nodesPositions[i + 1ull] << ' '
    //        << shell->nodesPositions[i + 2ull] << '\n';
    //polycr_mesh << '\n';
    //
    //polycr_mesh << shell->nFacets << '\n';
    //for (size_t i = 0ull; i < 3ull * shell->nFacets; i += 3ull)
    //{
    //    polycr_mesh
    //        << shell->facets[i] << ' '
    //        << shell->facets[i + 1ull] << ' '
    //        << shell->facets[i + 2ull] << '\n';
    //}
    //polycr_mesh << '\n';
    //
    //polycr_mesh << shell->nCryses << '\n';
    //for (size_t i = 0ull; i < shell->nCryses; i++)
    //    polycr_mesh << shell->nCrysesFacets[i] << ' ';
    //polycr_mesh << '\n';
    //
    //for (size_t i = 0ull; i < shell->nCryses; i++)
    //{
    //    for (int j = 0; j < 12; j++)
    //        polycr_mesh << shell->cryses[12ull * i + j] << ' ';
    //    polycr_mesh << '\n';
    //}
    
    std::cout << "\nInitializing polycrystal data...";
    Polycrystal3 polycr(shell);
    std::cout << " done.\n";
    delete shell;

    std::cout << "Generating mesh...";
    PolycrMesh* mesh = polycr.generateMesh(preferredEdgeLength);
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polycr.outputData("polycr_729_crysts.k", K);
    polycr.outputData("polycr_729_crysts.obj", OBJ);
    std::cout << " done.\n";

    #else

    std::cout << "\nInitializing polycrystal data...";
    Polycrystal3 polycr("input_cube_test.txt");
    std::cout << " done.\n";

    std::cout << "Generating mesh...";
    PolycrMesh* mesh = polycr.generateMesh(preferredEdgeLength); // 0.45 is quite good
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polycr.outputData();
    std::cout << " done.\n";

    #endif

    delete mesh;
    return 0;
}