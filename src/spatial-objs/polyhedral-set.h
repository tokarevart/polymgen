// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <memory>
// TODO: rewrite classes in STL style using templates and generalized algorithms working with that classes
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
#include "spatial-objs/genparams.h"
#include "spatial-objs/filetype.h"

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
        std::size_t nPolyhs   = 0;
        std::size_t nElems    = 0;
        real_t prefLen        = static_cast<real_t>(-1.0);
        double shellTrTime         = -1.0;
        double volumeExhTime       = -1.0;
        double meshFileWritingTime = -1.0;

        void write( std::ostream& stream ) const;
    };

    void tetrahedralize( real_t preferredLength, genparams::Polyhedron genParams = genparams::Polyhedron() );
    void smoothMesh( std::size_t nItersVolume = 1, std::size_t nItersShell = 1 );
    void shellDelaunayPostP();

    PolyMesh mesh() const;

    Log log();
    std::string generateLogFileName() const;

    // TODO: make mesh output with PolyMesh method instead
    void output( FileType filetype = FileType::WavefrontObj, std::string_view filename = "_AUTO_", unsigned polyhedralSetId = 1u );

    PolyhedralSet();
    PolyhedralSet( const psg::PolyShell& polyStruct );
    ~PolyhedralSet();


private:
    Log m_log;
    bool m_isLogged = false;

    real_t m_prefLen;

    std::vector<Polyhedron*> m_polyhedrons;

    // TODO: use Shell class instead
    std::vector<shell::Face*> m_shellFaces;
    std::vector<shell::Edge*> m_shellEdges;
    std::vector<shell::Vert*> m_shellVerts;

    shell::Edge* findShellEdge( const shell::Vert* v0, const shell::Vert* v1 ) const;

    void triangulateShell( genparams::Shell genParams );

    void outputObj( std::string_view filename ) const;
    void outputLSDynaKeyword_PART( std::ofstream& file ) const;
    void outputLSDynaKeyword_NODE( std::ofstream& file ) const;
    void outputLSDynaKeyword_ELEMENT_SOLID( std::ofstream& file, unsigned PolyhedralSetId = 1u ) const;
    void outputLSDynaKeyword( const std::string& filename, unsigned PolyhedralSetId = 1u ) const;

    std::string generateOutputFilename( FileType filetype, std::string_view filename ) const;

    void assign( std::string_view polyStructFileName ); // NOTE: deprecated
    void assign( const psg::PolyShell& polyStruct );
};

} // namespace pmg
