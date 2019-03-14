#include <stddef.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "spatial-objs/polycrystal.h"
#include "polygen/polygen.h"


int main()
{
    double preferredEdgeLength = 0.30;
//    std::cout << "Enter preferred shell::Edge length: ";
//    std::cin >> preferredEdgeLength;
    
    std::cout << "Generating polycrystal...";
    size_t n = 1;
    polygen::PolyStruct polystr = polygen::generateCuboidsPolycrystal(n, n, n);
    std::cout << " done.\n";

    std::cout << "Initializing polycrystal data...";
    pmg::Polycrystal polycr(polystr);
    std::cout << " done.\n";
    polystr.clear();

    std::cout << "Generating mesh...";
    polycr.generateMesh(preferredEdgeLength);
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polycr.output(pmg::Polycrystal::FileType::LsDynaKeyword);
    polycr.output(pmg::Polycrystal::FileType::Obj);
    std::cout << " done.\n";

    return 0;
}
