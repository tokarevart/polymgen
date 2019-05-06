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
#include "pmg-settings.h"
#include "file-type.h"

#include "definitions.h"


namespace pmg {

class PolyhedralSet
{
public:
    struct Log
    {
        real_t minQuality     = static_cast<real_t>(-1.0);
        real_t avQuality      = static_cast<real_t>(-1.0);
        real_t minMeshAbsGrad = static_cast<real_t>(-1.0);
        real_t avMeshAbsGrad  = static_cast<real_t>(-1.0);
        size_t nPolyhs        = 0;
        size_t nElems         = 0;
        real_t prefLen        = static_cast<real_t>(-1.0);
        double shellTrTime         = -1.0;
        double volumeExhTime       = -1.0;
        double meshFileWritingTime = -1.0;

        void write( std::ostream& stream ) const;
    };

    void tetrahedralize( real_t preferredLength, gensettings::Polyhedron genSettings = gensettings::Polyhedron() );
    void smoothMesh( size_t nItersVolume = 1, size_t nItersShell = 1 );
    void shellDelaunayPostP();
    const PolyMesh* structurizeMesh();
    const PolyMesh* getLastMesh();

    Log log();
    std::string generateLogFileName() const;

    void input( std::string_view polyStructFileName );
    void input( const psg::PolyShell& polyStruct );
    void output( FileType filetype = FileType::WavefrontObj, std::string_view filename = "_AUTO_", unsigned polyhedralSetId = 1u );

    PolyhedralSet();
    PolyhedralSet( std::string_view polyStructFileName );
    PolyhedralSet( const psg::PolyShell& polyStruct );
    ~PolyhedralSet();


private:
    Log m_log;
    bool m_isLogged = false;

    real_t m_prefLen;
    PolyMesh* m_lastMesh = nullptr;

    std::vector<Polyhedron*> m_polyhedrons;

    std::vector<shell::Face*> m_shellFaces;
    std::vector<shell::Edge*> m_shellEdges;
    std::vector<shell::Vert*> m_shellVerts;

    shell::Edge* findShellEdge( const shell::Vert* v0, const shell::Vert* v1 ) const;

    void triangulateShell( gensettings::Shell genSettings );

    void outputObj( std::string_view filename ) const;
    void outputLSDynaKeyword_PART( std::ofstream& file ) const;
    void outputLSDynaKeyword_NODE( std::ofstream& file ) const;
    void outputLSDynaKeyword_ELEMENT_SOLID( std::ofstream& file, unsigned PolyhedralSetId = 1u ) const;
    void outputLSDynaKeyword( const std::string& filename, unsigned PolyhedralSetId = 1u ) const;

    std::string generateOutputFilename( FileType filetype, std::string_view filename ) const;
};

} // namespace pmg
