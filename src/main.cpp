// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include <stddef.h>
#include <iostream>
#include <fstream>
#include <memory>
#include "spatial-objs/polyhedral-set.h"
#include "polysgen/polysgen.h"


int main()
{
    real_t preferredEdgeLength = static_cast<real_t>(0.45);
    
    std::cout << "Generating PolyShell...";
    size_t n = 4;
    psg::PolyShell polys = psg::generateCuboids(n, n, n);
    std::cout << " done.\n";

    std::cout << "Initializing PolyhedralSet...";
    pmg::PolyhedralSet polyhedr(polys);
    std::cout << " done.\n";
    polys.clear();

    std::cout << "Generating mesh...";
    polyhedr.generateMesh(preferredEdgeLength);
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polyhedr.output(pmg::PolyhedralSet::FileType::LsDynaKeyword);
    polyhedr.output(pmg::PolyhedralSet::FileType::Obj);
    std::cout << " done.\n";

    return 0;
}
