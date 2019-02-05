#include <stddef.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "Polycrystal3.h"
#include "polygen/PolyGen.h"


int main()
{
    double preferredEdgeLength = 0.45;
//    std::cout << "Enter preferred edge length: ";
//    std::cin >> preferredEdgeLength;
    
    std::cout << "Generating polycrystal...";
    size_t n = 4;
    PolyStruct polystr = polygen::generateCuboidsPolycrystal(n, n, n);
    std::cout << " done.\n";

    std::cout << "\nInitializing polycrystal data...";
    Polycrystal3 polycr(polystr);
    std::cout << " done.\n";
    polystr.clear();

    std::cout << "Generating mesh...";
    polycr.generateMesh(preferredEdgeLength);
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polycr.outputData(Polycrystal3::FileType::LS_DYNA_KEYWORD);
//    polycr.outputData(Polycrystal3::OBJ);
    std::cout << " done.\n";

    return 0;
}
