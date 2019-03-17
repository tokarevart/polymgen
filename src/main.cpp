// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Contacts: <tokarev28.art@gmail.com>
// Licensed under the MIT License.

#include <stddef.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "spatial-objs/polycrystal.h"
#include "polygen/polygen.h"


// Note:
//   Try to remove EPS calculation in spatalgs::distancePointToSegment
//   and spatalgs::closestSegmentPointToPoint.
//   Try using my new compare functions.

int main()
{
    double preferredEdgeLength = 0.45;
//    std::cout << "Enter preferred shell::Edge length: ";
//    std::cin >> preferredEdgeLength;
    
    std::cout << "Generating polycrystal...";
    size_t n = 3;
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
