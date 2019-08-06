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


std::string rabbit() {
    std::array<std::string, 4> lines =
    {
        "(\\_/)",
        "(._.)",
        "/   \\",
        "-----"
    };
    std::string res;
    for (auto&& line : lines)
        res += line + '\n';
    return res;
}
#include "src/helpers/spatial/algs.h"
int main() {
    std::cout << rabbit() << std::endl;

    real_t preferred_length = static_cast<real_t>(0.17);

    std::cout << "Generating PolyShell...";
    std::size_t n = 4;
    psg::PolyShell shell = psg::cuboids(n, n, n);
    std::cout << std::string(7, ' ') + "done." << std::endl;

    std::cout << "Initializing PolyhedralSet...";
    pmg::PolyhedralSet polyhset(shell);
    std::cout << " done." << std::endl;
    shell.clear();

    std::cout << "Generating mesh...";
    polyhset.tetrahedralize(preferred_length);
    std::cout << std::string(12, ' ') + "done." << std::endl;

    //std::cout << std::endl << polyhset.log_file_name();

    namespace fs = std::filesystem;
    std::cout << "Outputting data to file...";

    fs::create_directories(fs::current_path() / "out" / "meshes");
    fs::current_path(fs::current_path() / "out" / "meshes");

    polyhset.output(pmg::filetype::lsdyna_keyword);
    polyhset.output(pmg::filetype::wavefront_obj);
    std::cout << std::string(4, ' ') + "done." << std::endl << std::endl;

    fs::current_path(fs::current_path().parent_path());
    fs::create_directory(fs::current_path() / "logs");
    fs::current_path(fs::current_path() / "logs");

    std::ofstream log_file(polyhset.log_file_name());
    polyhset.log().write(log_file);
    std::cout << "Log:" << std::endl;
    polyhset.log().write(std::cout);
    log_file.close();

    return 0;
}
