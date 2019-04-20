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
    real_t preferredLength = static_cast<real_t>(0.5);
    
    std::cout << "Generating PolyShell...";
    size_t n = 4;
    psg::PolyShell shell = psg::generateCuboids(n, n, n);
    std::cout << " done.\n";

    std::cout << "Initializing PolyhedralSet...";
    pmg::PolyhedralSet polyhset(shell);
    std::cout << " done.\n";
    shell.clear();

    std::cout << "Generating mesh...";
    polyhset.generateMesh(preferredLength);
    polyhset.optimizeMesh();
    std::cout << " done.\n";

    std::cout << "Outputting data to file...";
    polyhset.output(pmg::PolyhedralSet::FileType::LsDynaKeyword);
    polyhset.output(pmg::PolyhedralSet::FileType::WavefrontObj);
    std::cout << " done.\n";

    std::ofstream log_file(polyhset.generateLogFileName());
    polyhset.log().write(log_file);
    log_file.close();

    return 0;
}
