#include <iostream>
#include <fstream>
#include "Polycrystal3.h"
#include "PolyGen.h"

int main(int argc, char* argv[])
{
    double preferredEdgeLength = 0.45;
    //std::cout << "Enter preferred edge length: ";
    //std::cin >> preferredEdgeLength;

    #define GENERATE
    #ifdef GENERATE

    std::cout << "Generating polycrystal...";
    PolyStruct* polystr = polygen::generateCuboidsPolycrystal(6, 6, 6);
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
    Polycrystal3 polycr(polystr);
    std::cout << " done.\n";
    delete polystr;

    std::cout << "Generating mesh...";
    polycr.generateMeshNoStructGen(preferredEdgeLength);
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polycr.outputData("polycr_216_crysts.kw", LS_DYNA_KEYWORD);
    //polycr.outputData("polycr_216_crysts.obj", OBJ);
    std::cout << " done.\n";

    #else

    std::cout << "\nInitializing polycrystal data...";
    Polycrystal3 polycr("input_cube_test.txt");
    std::cout << " done.\n";

    std::cout << "Generating mesh...";
    PolyMesh* mesh = polycr.generateMesh(preferredEdgeLength); // 0.45 is quite good
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polycr.outputData();
    std::cout << " done.\n";

    #endif

    return 0;
}