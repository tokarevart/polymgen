// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include <cstddef>
#include <iostream>
#include <fstream>
#include <memory>
#include <string>
#include <filesystem>
#include "src/core/polyhedral-set.h"
#include "src/polysgen/polysgen.h"

int main()
{
    real_t preferredLength = static_cast<real_t>(0.5);
    
    std::cout << "Generating PolyShell...";
    std::size_t n = 4;
    psg::PolyShell shell = psg::generateCuboids(n, n, n);
    std::cout << std::string(7, ' ') + "done." << std::endl;

    std::cout << "Initializing PolyhedralSet...";
    pmg::PolyhedralSet polyhset(shell);
    std::cout << " done." << std::endl;
    shell.clear();

    std::cout << "Generating mesh...";
    polyhset.tetrahedralize(preferredLength);
    std::cout << std::string(12, ' ') + "done." << std::endl;

    namespace fs = std::filesystem;
    std::cout << "Outputting data to file...";

    fs::create_directories(fs::current_path() / "out" / "meshes");
    fs::current_path(fs::current_path() / "out" / "meshes");

    polyhset.output(pmg::FileType::LsDynaKeyword);
    polyhset.output(pmg::FileType::WavefrontObj);
    std::cout << std::string(4, ' ') + "done." << std::endl << std::endl;

    fs::current_path(fs::current_path().parent_path());
    fs::create_directory(fs::current_path() / "logs");
    fs::current_path(fs::current_path() / "logs");

    std::ofstream log_file(polyhset.generateLogFileName());
    polyhset.log().write(log_file);
    std::cout << "Log:" << std::endl;
    polyhset.log().write(std::cout);
    log_file.close();

    return 0;
}
