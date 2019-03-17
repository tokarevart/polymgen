// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include "spatial-objs/crystallite.h"
#include "spatial-objs/shell/shell-facet.h"
#include "spatial-objs/shell/shell-edge.h"
#include "spatial-objs/shell/shell-vertex.h"
#include "spatial-objs/tetr.h"
#include "spatial-objs/facet.h"
#include "spatial-objs/edge.h"
#include "spatial-objs/vertex.h"
#include "data-structures/polymesh.h"
#include "data-structures/polystruct.h"
#include "helpers/logger.h"

#include "definitions.h"


namespace pmg {

class Polycrystal
{
public:
    enum class FileType
    {
        Obj,
        LsDynaKeyword
    };

    void generateMesh(double preferredLength, std::string_view logFileName = "_AUTO_");
    const PolyMesh* structurizeMesh();
    const PolyMesh* getLastMesh();

    void input(std::string_view polyStructFileName);
    void input(const polygen::PolyStruct& polyStruct);
    void output(FileType filetype = FileType::Obj, std::string_view filename = "_AUTO_", unsigned polycrystalId = 1u) const;

    Polycrystal();
    Polycrystal(std::string_view polyStructFileName);
    Polycrystal(const polygen::PolyStruct& polyStruct);
    ~Polycrystal();


private:
    double m_preferredLength;

    PolyMesh* m_lastMesh = nullptr;
    std::unique_ptr<tva::Logger> m_lastLogger;

    std::vector<Crystallite*> m_crystallites;

    std::vector<shell::Facet*>  m_shellFacets;
    std::vector<shell::Edge*>   m_shellEdges;
    std::vector<shell::Vertex*> m_shellVerts;

    shell::Edge* findShellEdge(const shell::Vertex* v0, const shell::Vertex* v1) const;

    void triangulateShell();

    void outputObj(std::string_view filename) const;
    void outputLSDynaKeyword_PART(std::ofstream& file) const;
    void outputLSDynaKeyword_NODE(std::ofstream& file) const;
    void outputLSDynaKeyword_ELEMENT_SOLID(std::ofstream& file, unsigned polycrystalId = 1u) const;
    void outputLSDynaKeyword(const std::string& filename, unsigned polycrystalId = 1u) const;

    std::string generateLogFileName(std::string_view logFileName) const;
    std::string generateOutputFilename(FileType filetype, std::string_view filename) const;
};

} // namespace pmg
