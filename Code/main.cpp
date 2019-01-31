#include <iostream>
#include <fstream>
#include "Polycrystal3.h"
#include "PolyGen/PolyGen.h"


int main(int argc, char* argv[])
{
    double preferredEdgeLength = 0.45;
    //std::cout << "Enter preferred edge length: ";
    //std::cin >> preferredEdgeLength;
    
    std::cout << "Generating polycrystal...";
    int n = 4;
    PolyStruct* polystr = polygen::generateCuboidsPolycrystal(n, n, n);
    std::cout << " done.\n";
        
    std::cout << "\nInitializing polycrystal data...";
    Polycrystal3 polycr(polystr);
    std::cout << " done.\n";
    delete polystr;

    std::cout << "Generating mesh...";
    polycr.generateMesh(preferredEdgeLength);
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polycr.outputData(Polycrystal3::LS_DYNA_KEYWORD);
    //polycr.outputData(Polycrystal3::OBJ);
    std::cout << " done.\n";

    return 0;
}