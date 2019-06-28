// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#pragma once
#include <cstdint>
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

// make STL naming in all the library
namespace pmg {

// refactor class
class PolyhedralSet {
public:
    struct Log {
        real_t min_quality = static_cast<real_t>(-1.0);
        real_t av_quality = static_cast<real_t>(-1.0);
        real_t min_mesh_abs_grad = static_cast<real_t>(-1.0);
        real_t av_mesh_abs_grad = static_cast<real_t>(-1.0);
        std::size_t n_polyhs = 0;
        std::size_t n_elems = 0;
        real_t pref_len = static_cast<real_t>(-1.0);
        double shell_tr_time = -1.0;
        double volume_exh_time = -1.0;
        double mesh_file_writing_time = -1.0;

        void write(std::ostream& stream) const;
    };

    void tetrahedralize(real_t preferred_length, genparams::Polyhedron gen_params = genparams::Polyhedron());
    void smooth_mesh(std::size_t nItersVolume = 1, std::size_t nItersShell = 1);
    void shell_delaunay_postp();

    PolyMesh mesh() const;

    Log log();
    std::string log_file_name() const;

    // TODO: make mesh output with PolyMesh method instead
    // TODO: make method taking meshes of all the polyhedrons and remeshing them into one
    void output(FileType filetype = FileType::WavefrontObj, std::string_view filename = "_AUTO_", unsigned polyhedral_set_id = 1u);

    PolyhedralSet() = delete;
    PolyhedralSet(const psg::PolyShell& poly_struct);
    ~PolyhedralSet();


private:
    Log m_log;
    bool m_is_logged = false;

    real_t m_prefLen;

    std::vector<Polyhedron*> m_polyhedrons;

    // TODO: use Shell class instead
    std::vector<shell::Face*> m_shell_faces;
    std::vector<shell::Edge*> m_shell_edges;
    std::vector<shell::Vert*> m_shell_verts;

    shell::Edge* find_shell_edge(const shell::Vert* v0, const shell::Vert* v1) const;

    void triangulate_shell(genparams::Shell gen_params);

    void output_obj(std::string_view filename) const;
    void output_lsdynakw_part_section(std::ofstream& file) const;
    void output_lsdynakw_node_section(std::ofstream& file) const;
    void output_lsdynakw_element_solid_section(std::ofstream& file, std::uint8_t PolyhedralSetId = 1u) const;
    void output_lsdynakw(const std::string& filename, std::uint8_t PolyhedralSetId = 1u) const;

    std::string output_filename(FileType filetype, std::string_view filename) const;

    void set_poly_shell(std::string_view poly_struct_filename); // NOTE: deprecated
    void set_poly_shell(const psg::PolyShell& poly_struct);
};

} // namespace pmg
