// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <string>
#include <fstream>
#include <list>
#include <vector>
#include <memory>
#include "polyhedron.h"
#include "shell/face.h"
#include "shell/edge.h"
#include "shell/vert.h"
#include "tetr.h"
#include "face.h"
#include "edge.h"
#include "vert.h"
#include "../data-structs/polymesh.h"
#include "../data-structs/polyshell.h"
#include "../helpers/logger.h"
#include "../real-type.h"
#include "genparams.h"
#include "filetype.h"

#include "../definitions.h"


namespace pmg {
// refactor class
class PolyhedralSet // TODO: rename to PolyMesher<3> or PolyMesher3
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

    PolyhedralSet() = delete;
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

    void setPolyShell( std::string_view polyStructFileName ); // NOTE: deprecated
    void setPolyShell( const psg::PolyShell& polyStruct );
};

} // namespace pmg
