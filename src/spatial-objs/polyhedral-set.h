// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include "spatial-objs/polyhedron.h"
#include "spatial-objs/shell/shell-face.h"
#include "spatial-objs/shell/shell-edge.h"
#include "spatial-objs/shell/shell-vertex.h"
#include "spatial-objs/tetr.h"
#include "spatial-objs/face.h"
#include "spatial-objs/edge.h"
#include "spatial-objs/vert.h"
#include "data-structures/polymesh.h"
#include "data-structures/polyshell.h"
#include "helpers/logger.h"
#include "real-type.h"

#include "definitions.h"


namespace pmg {

class PolyhedralSet
{
public:
    enum class FileType
    {
        WavefrontObj,
        LsDynaKeyword
    };

    void generateMesh( real_t preferredLength, std::string_view logFileName = "_AUTO_" );
    const PolyMesh* structurizeMesh();
    const PolyMesh* getLastMesh();

    void input( std::string_view polyStructFileName );
    void input( const psg::PolyShell& polyStruct );
    void output( FileType filetype = FileType::WavefrontObj, std::string_view filename = "_AUTO_", unsigned PolyhedralSetId = 1u ) const;

    PolyhedralSet();
    PolyhedralSet( std::string_view polyStructFileName );
    PolyhedralSet( const psg::PolyShell& polyStruct );
    ~PolyhedralSet();


private:
    real_t m_preferredLength;

    PolyMesh* m_lastMesh = nullptr;
    std::unique_ptr<Logger> m_lastLogger;

    std::vector<Polyhedron*> m_polyhedrons;

    std::vector<shell::Face*> m_shellFaces;
    std::vector<shell::Edge*> m_shellEdges;
    std::vector<shell::Vert*> m_shellVerts;

    shell::Edge* findShellEdge( const shell::Vert* v0, const shell::Vert* v1 ) const;

    void triangulateShell();

    void outputObj( std::string_view filename ) const;
    void outputLSDynaKeyword_PART( std::ofstream& file ) const;
    void outputLSDynaKeyword_NODE( std::ofstream& file ) const;
    void outputLSDynaKeyword_ELEMENT_SOLID( std::ofstream& file, unsigned PolyhedralSetId = 1u ) const;
    void outputLSDynaKeyword( const std::string& filename, unsigned PolyhedralSetId = 1u ) const;

    std::string generateLogFileName( std::string_view logFileName ) const;
    std::string generateOutputFilename( FileType filetype, std::string_view filename ) const;
};

} // namespace pmg
