#pragma once
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include "Crystallite3.h"
#include "Simplex3.h"
#include "Facet3.h"
#include "ShellFacet3.h"
#include "Edge3.h"
#include "ShellEdge3.h"
#include "Vertex3.h"
#include "ShellVertex3.h"
#include "data-structures/PolyMesh.h"
#include "data-structures/PolyStruct.h"
#include "helpers/Logger.h"

class Crystallite3;
class Simplex3;
class FrontFacet3;
class ShellFacet3;
class Facet3;
class FrontEdge3;
class ShellEdge3;
class Edge3;
class Vertex3;
class ShellVertex3;


class Polycrystal3
{
public:
    enum FileType
    {
        OBJ,
        LS_DYNA_KEYWORD
    };

    void generateMesh(const double preferredLength, const std::string& logFileName = "_AUTO_");
    void generateMesh(const std::string& polyStructFileName, const double preferredLength, const std::string& logFileName = "_AUTO_");
    void generateMesh(const std::unique_ptr<PolyStruct>& polyStruct, const double preferredLength, const std::string& logFileName = "_AUTO_");
    const PolyMesh* structurizeMesh();
    const PolyMesh* getLastMesh();

    void inputData(const std::string& polyStructFileName);
    void inputData(const std::unique_ptr<PolyStruct>& polyStruct);
    void outputData(FileType filetype = OBJ, const std::string& filename = "_AUTO_", unsigned polycrystalId = 1u) const;

    Polycrystal3();
    Polycrystal3(const std::string& polyStructFileName);
    Polycrystal3(const std::unique_ptr<PolyStruct>& polyStruct);
    ~Polycrystal3();


private:
    double m_preferredLength;

    PolyMesh* m_lastMesh = nullptr;
    std::unique_ptr<tva::Logger> m_lastLogger;

    std::vector<Crystallite3*> m_crystallites;

    std::vector<ShellFacet3*> m_shellFacets;
    std::vector<ShellEdge3*> m_shellEdges;
    std::vector<ShellVertex3*> m_shellVertexes;

    std::list<Facet3*> m_startFrontFacets;
    std::list<Edge3*> m_startFrontEdges;
    std::list<Vertex3*> m_startFrontVertexes;

    // Later replace ShellVertex3::findAttachedVertex with that function
    Vertex3* findAttachedVertex(const ShellVertex3* shellVertex);
    ShellEdge3* findShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const;
    Edge3* findStartFrontEdge(const Vertex3* v0, const Vertex3* v1) const;

    void triangulateShell();
    void setLinksWithShell();
    void startFrontDelaunayPostprocessing();

    void outputDataObj(const std::string& filename) const;
    void outputDataLSDynaKeyword_PART(std::ofstream& file) const;
    void outputDataLSDynaKeyword_NODE(std::ofstream& file) const;
    void outputDataLSDynaKeyword_ELEMENT_SOLID(std::ofstream& file, unsigned polycrystalId = 1u) const;
    void outputDataLSDynaKeyword(const std::string& filename, unsigned polycrystalId = 1u) const;

    std::string generateMesh_generateLogFileName(const std::string& logFileName) const;
    std::string outputData_generateFilename(FileType filetype, const std::string& filename) const;
};
