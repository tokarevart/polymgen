// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include <stddef.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "spatial-objs/polyhedral-set.h"
#include "polygen/polygen.h"


int main()
{
    real_t preferredEdgeLength = static_cast<real_t>(0.16);
//    std::cout << "Enter preferred shell::Edge length: ";
//    std::cin >> preferredEdgeLength;
    
    std::cout << "Generating PolyhedralSet...";
    size_t n = 1;
    polygen::PolyStruct polystr = polygen::generateCuboidsPolyhedralSet(n, n, n);
    std::cout << " done.\n";

    std::cout << "Initializing PolyhedralSet data...";
    pmg::PolyhedralSet polycr(polystr);
    std::cout << " done.\n";
    polystr.clear();

    std::cout << "Generating mesh...";
    polycr.generateMesh(preferredEdgeLength);
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polycr.output(pmg::PolyhedralSet::FileType::LsDynaKeyword);
    polycr.output(pmg::PolyhedralSet::FileType::Obj);
    std::cout << " done.\n";

    return 0;
}
