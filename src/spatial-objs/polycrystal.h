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

    void generateMesh(double preferredLength, const std::string& logFileName = "_AUTO_");
    const PolyMesh* structurizeMesh();
    const PolyMesh* getLastMesh();

    void input(const std::string& polyStructFileName);
    void input(const polygen::PolyStruct& polyStruct);
    void output(FileType filetype = FileType::Obj, const std::string& filename = "_AUTO_", unsigned polycrystalId = 1u) const;

    Polycrystal();
    Polycrystal(const std::string& polyStructFileName);
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

    void outputObj(const std::string& filename) const;
    void outputLSDynaKeyword_PART(std::ofstream& file) const;
    void outputLSDynaKeyword_NODE(std::ofstream& file) const;
    void outputLSDynaKeyword_ELEMENT_SOLID(std::ofstream& file, unsigned polycrystalId = 1u) const;
    void outputLSDynaKeyword(const std::string& filename, unsigned polycrystalId = 1u) const;

    std::string generateMesh_generateLogFileName(const std::string& logFileName) const;
    std::string output_generateFilename(FileType filetype, const std::string& filename) const;
};

} // namespace pmg
