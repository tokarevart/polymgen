// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "polyhedral-set.h"
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <chrono>
#include "../helpers/spatial/algs.h"

using namespace pmg;
using vec3 = spt::vec<3, real_t>;


#define DEBUG


shell::Edge* PolyhedralSet::find_shell_edge(const shell::Vert* v0, const shell::Vert* v1) const {
    for (auto& s_edge : m_shell_edges)
        if ((s_edge->verts[0] == v0
             && s_edge->verts[1] == v1)
            || (s_edge->verts[1] == v0
                && s_edge->verts[0] == v1))
            return s_edge;

    return nullptr;
}


void PolyhedralSet::triangulate_shell(genparams::Shell gen_params) {
    for (auto& svert : m_shell_verts)
        svert->attached_vert = new pmg::Vert(svert->pos());

    for (auto& sedge : m_shell_edges)
        sedge->segmentize(m_pref_len);

    std::size_t n_shell_faces = m_shell_faces.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_shell_faces; i++)
        m_shell_faces[i]->triangulate(m_pref_len, gen_params);
}


void PolyhedralSet::output_obj(std::string_view filename) const {
    std::ofstream file(filename.data());

    std::size_t index = 1;
    for (auto& svert : m_shell_verts) {
        file << "v " << svert->pos()[0] << ' ' << svert->pos()[1] << ' ' << svert->pos()[2] << '\n';
        svert->attached_vert->global_idx = index++;
    }
    for (auto& sedge : m_shell_edges) {
        for (auto& vert : sedge->inner_verts()) {
            file << "v " << vert->pos()[0] << ' ' << vert->pos()[1] << ' ' << vert->pos()[2] << '\n';
            vert->global_idx = index++;
        }
    }
    for (auto& sface : m_shell_faces) {
        for (auto& vert : sface->inner_verts()) {
            file << "v " << vert->pos()[0] << ' ' << vert->pos()[1] << ' ' << vert->pos()[2] << '\n';
            vert->global_idx = index++;
        }
    }

    for (auto& polyhedr : m_polyhedrons) {
        for (auto& vert : polyhedr->inner_verts()) {
            file << "v " << vert->pos()[0] << ' ' << vert->pos()[1] << ' ' << vert->pos()[2] << '\n';
            vert->global_idx = index++;
        }
    }

    #ifdef DEBUG
    for (auto& polyhedr : m_polyhedrons) {
        for (auto& f_face : polyhedr->front_faces()) {
            auto face = f_face->x;

            std::vector<std::size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->global_idx) == gl_nums.end())
                        gl_nums.push_back(vert->global_idx);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
    #else
    for (auto& sface : m_shell_faces) {
        for (auto& face : sface->inner_faces()) {
            std::vector<std::size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->global_idx) == gl_nums.end())
                        gl_nums.push_back(vert->global_idx);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }

    for (auto& polyhedr : m_polyhedrons) {
        for (auto& face : polyhedr->inner_faces()) {
            std::vector<std::size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->global_idx) == gl_nums.end())
                        gl_nums.push_back(vert->global_idx);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
    #endif

    file.close();
}


void PolyhedralSet::output_lsdynakw_part_section(std::ofstream& file) const {
    for (std::size_t id = 1; id < m_polyhedrons.size() + 1; id++) {
        file << "*PART\n"
            "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n";
        file << std::setw(10) << id;
        file << "         0";
        file << std::setw(10) << id;
        file << "         0         0         0         0         0\n";
    }
}


void PolyhedralSet::output_lsdynakw_node_section(std::ofstream& file) const {
    file << "*NODE\n"
        "$#   nid               x               y               z      tc      rc\n";
    std::size_t index = 1;
    for (auto& svert : m_shell_verts) {
        svert->attached_vert->global_idx = index++;
        file << std::setw(8) << svert->attached_vert->global_idx;
        file << std::setw(16) << svert->pos()[0];
        file << std::setw(16) << svert->pos()[1];
        file << std::setw(16) << svert->pos()[2];
        file << "       0       0\n";
    }
    for (auto& sedge : m_shell_edges) {
        for (auto& vert : sedge->inner_verts()) {
            vert->global_idx = index++;
            file << std::setw(8) << vert->global_idx;
            file << std::setw(16) << vert->pos()[0];
            file << std::setw(16) << vert->pos()[1];
            file << std::setw(16) << vert->pos()[2];
            file << "       0       0\n";
        }
    }
    for (auto& sface : m_shell_faces) {
        for (auto& vert : sface->inner_verts()) {
            vert->global_idx = index++;
            file << std::setw(8) << vert->global_idx;
            file << std::setw(16) << vert->pos()[0];
            file << std::setw(16) << vert->pos()[1];
            file << std::setw(16) << vert->pos()[2];
            file << "       0       0\n";
        }
    }

    for (auto& polyhedr : m_polyhedrons) {
        for (auto& vert : polyhedr->inner_verts()) {
            vert->global_idx = index++;
            file << std::setw(8) << vert->global_idx;
            file << std::setw(16) << vert->pos()[0];
            file << std::setw(16) << vert->pos()[1];
            file << std::setw(16) << vert->pos()[2];
            file << "       0       0\n";
        }
    }
}


void PolyhedralSet::output_lsdynakw_element_solid_section(std::ofstream& file, std::uint8_t PolyhedralSetId) const {
    file << "*ELEMENT_SOLID\n"
        "$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n";
    std::size_t pid = 1;
    std::size_t eid = 10000000ull * PolyhedralSetId + 1ull;
    for (auto& polyhedr : m_polyhedrons) {
        for (auto& tetr : polyhedr->inner_tetrs()) {
            file << std::setw(8) << eid++;
            file << std::setw(8) << pid;
            file << std::setw(8) << tetr->verts[0]->global_idx;
            vec3 v0 = tetr->verts[1]->pos() - tetr->verts[0]->pos();
            vec3 v1 = tetr->verts[2]->pos() - tetr->verts[0]->pos();
            vec3 v2 = tetr->verts[3]->pos() - tetr->verts[0]->pos();
            if (spt::dot(v2, spt::cross(v0, v1)) > static_cast<real_t>(0.0)) {
                file << std::setw(8) << tetr->verts[1]->global_idx;
                file << std::setw(8) << tetr->verts[2]->global_idx;
            } else {
                file << std::setw(8) << tetr->verts[2]->global_idx;
                file << std::setw(8) << tetr->verts[1]->global_idx;
            }
            file << std::setw(8) << tetr->verts[3]->global_idx;
            file << "       0       0       0       0\n";
        }
        pid++;
    }
}


void PolyhedralSet::output_lsdynakw(const std::string& filename, std::uint8_t PolyhedralSetId) const {
    std::ofstream file(filename);
    file.setf(std::ios::right | std::ios::fixed);

    output_lsdynakw_part_section(file);
    output_lsdynakw_node_section(file);
    output_lsdynakw_element_solid_section(file, PolyhedralSetId);
    file << "*END";

    file.close();
}


std::string PolyhedralSet::log_file_name() const {
    std::stringstream ss;
    ss << "pmg_log_" << m_polyhedrons.size() << "_nph_";
    std::size_t av_nfe = 0;
    for (auto& polyhedr : m_polyhedrons)
        av_nfe += polyhedr->inner_tetrs().size();
    av_nfe /= m_polyhedrons.size();
    ss << av_nfe << "_phfe.log";

    return ss.str();
}


void PolyhedralSet::tetrahedralize(real_t preferred_length, genparams::Polyhedron gen_params) {
    m_pref_len = preferred_length;

    using namespace std::chrono;
    auto shell_triang_start = steady_clock::now();
    try {
        triangulate_shell(gen_params.shell);
    } catch (std::logic_error& error) {
        output(filetype::wavefront_obj, "debug.obj");
        throw error;
    }
    auto shell_triang_stop = steady_clock::now();
    duration<double> shell_triang_elapsed = shell_triang_stop - shell_triang_start;

    auto volume_exh_start = steady_clock::now();
    std::size_t n_polyhs = m_polyhedrons.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_polyhs; i++) {
        try {
            m_polyhedrons[i]->tetrahedralize(m_pref_len, gen_params.volume);
        } catch (std::logic_error& error) {
            output(filetype::wavefront_obj, "debug.obj");
            throw error;
        }
    }
    auto volume_exh_stop = steady_clock::now();
    duration<double> volume_exh_elapsed = volume_exh_stop - volume_exh_start;

    m_log.n_polyhs = m_polyhedrons.size();
    m_log.pref_len = m_pref_len;
    m_log.shell_tr_time = shell_triang_elapsed.count();
    m_log.volume_exh_time = volume_exh_elapsed.count();
    m_is_logged = false;
}


void PolyhedralSet::smooth_mesh(std::size_t n_iters_volume, std::size_t n_iters_shell) {
    std::size_t n_sfaces = m_shell_faces.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_sfaces; i++)
        m_shell_faces[i]->smooth_mesh(n_iters_shell);

    std::size_t n_polyhs = m_polyhedrons.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_polyhs; i++)
        m_polyhedrons[i]->smooth_mesh(n_iters_volume);

    m_is_logged = false;
}


void PolyhedralSet::shell_delaunay_postp() {
    std::size_t n_sfaces = m_shell_faces.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_sfaces; i++)
        m_shell_faces[i]->delaunay_postp();

    m_is_logged = false;
}


PolyMesh PolyhedralSet::mesh() const {
    PolyMesh mesh;
    mesh.polyhs.assign(m_polyhedrons.size(), std::vector<PolyMesh::TetrIdx>());

    std::size_t n_verts = m_shell_verts.size();
    for (auto& sedge : m_shell_edges)
        n_verts += sedge->inner_verts().size();
    for (auto& sface : m_shell_faces)
        n_verts += sface->inner_verts().size();
    for (auto& polyh : m_polyhedrons)
        n_verts += static_cast<std::size_t>(std::count_if(
            polyh->inner_verts().begin(),
            polyh->inner_verts().end(),
            [](auto vert) { return static_cast<bool>(vert); }));

    mesh.verts.assign(n_verts, PolyMesh::Vert());

    std::size_t vert_idx = 0;
    for (auto& svert : m_shell_verts) {
        mesh.verts[vert_idx] = svert->pos().x;
        svert->attached_vert->global_idx = vert_idx;
        vert_idx++;
    }
    auto add_verts_to_last_mesh = [&mesh, &vert_idx](auto verts) mutable {
        for (auto& vert : verts) {
            mesh.verts[vert_idx] = vert->pos().x;
            vert->global_idx = vert_idx;
            vert_idx++;
        }
    };
    for (auto& sedge : m_shell_edges) add_verts_to_last_mesh(sedge->inner_verts());
    for (auto& sface : m_shell_faces) add_verts_to_last_mesh(sface->inner_verts());
    for (auto& polyh : m_polyhedrons) add_verts_to_last_mesh(polyh->inner_verts());

    std::size_t n_tetrs = 0;
    for (std::size_t i = 0; i < mesh.polyhs.size(); i++) {
        mesh.polyhs[i].assign(m_polyhedrons[i]->inner_tetrs().size(), 0);
        n_tetrs += mesh.polyhs[i].size();
    }

    mesh.tetrs.assign(n_tetrs, { 0, 0, 0, 0 });

    std::size_t tetr_idx = 0;
    for (auto& polyh : m_polyhedrons)
        for (auto& tetr : polyh->inner_tetrs()) {
            std::array<PolyMesh::VertIdx, 4> verts_idces;
            verts_idces = { tetr->verts[0]->global_idx,
                            tetr->verts[1]->global_idx,
                            tetr->verts[2]->global_idx,
                            tetr->verts[3]->global_idx };

            vec3 v0 = tetr->verts[1]->pos() - tetr->verts[0]->pos();
            vec3 v1 = tetr->verts[2]->pos() - tetr->verts[0]->pos();
            vec3 v2 = tetr->verts[3]->pos() - tetr->verts[0]->pos();
            if (spt::dot(v2, spt::cross(v0, v1)) < static_cast<real_t>(0.0))
                std::swap(verts_idces[0], verts_idces[1]);

            mesh.tetrs[tetr_idx++] = verts_idces;
        }

    return mesh;
}


PolyhedralSet::Log PolyhedralSet::log() {
    if (m_is_logged)
        return m_log;

    std::size_t n_elems = 0;
    real_t min_q = 1.0, av_q = 0.0;
    real_t min_g = 1.0, av_g = 0.0;
    for (auto& polyh : m_polyhedrons) {
        auto cur_min_av_q = polyh->analyze_mesh_quality();
        auto cur_min_av_g = polyh->analyze_mesh_abs_grad();
        min_q = std::min(cur_min_av_q.first, min_q);
        min_g = std::min(cur_min_av_g.first, min_g);
        av_q += cur_min_av_q.second;
        av_g += cur_min_av_g.second;
        n_elems += polyh->inner_tetrs().size();
    }
    av_q /= m_polyhedrons.size();
    av_g /= m_polyhedrons.size();

    m_log.min_quality = min_q;
    m_log.av_quality = av_q;
    m_log.min_mesh_abs_grad = min_g;
    m_log.av_mesh_abs_grad = av_g;
    m_log.n_elems = n_elems;

    return m_log;
}


void PolyhedralSet::set_poly_shell(std::string_view poly_struct_filename) {
    std::ifstream input(poly_struct_filename.data());

    std::size_t n_verts;
    input >> n_verts;

    for (std::size_t i = 0; i < n_verts; i++) {
        vec3 v;
        input >> v.x[0];
        input >> v.x[1];
        input >> v.x[2];

        m_shell_verts.push_back(new shell::Vert(v));
    }

    std::size_t faces_num;
    input >> faces_num;

    for (std::size_t i = 0; i < faces_num; i++) {
        std::array<std::size_t, 3> face_nodes_inds;
        input >> face_nodes_inds[0];
        input >> face_nodes_inds[1];
        input >> face_nodes_inds[2];

        if (!find_shell_edge(m_shell_verts[face_nodes_inds[0]], m_shell_verts[face_nodes_inds[1]])) {
            m_shell_edges.push_back(new shell::Edge(m_shell_verts[face_nodes_inds[0]], m_shell_verts[face_nodes_inds[1]]));
        }

        if (!find_shell_edge(m_shell_verts[face_nodes_inds[1]], m_shell_verts[face_nodes_inds[2]])) {
            m_shell_edges.push_back(new shell::Edge(m_shell_verts[face_nodes_inds[1]], m_shell_verts[face_nodes_inds[2]]));
        }

        if (!find_shell_edge(m_shell_verts[face_nodes_inds[2]], m_shell_verts[face_nodes_inds[0]])) {
            m_shell_edges.push_back(new shell::Edge(m_shell_verts[face_nodes_inds[2]], m_shell_verts[face_nodes_inds[0]]));
        }

        m_shell_faces.push_back(new shell::Face(
            find_shell_edge(m_shell_verts[face_nodes_inds[0]], m_shell_verts[face_nodes_inds[1]]),
            find_shell_edge(m_shell_verts[face_nodes_inds[1]], m_shell_verts[face_nodes_inds[2]]),
            find_shell_edge(m_shell_verts[face_nodes_inds[2]], m_shell_verts[face_nodes_inds[0]])));
    }

    std::size_t n_polyhedrons;
    input >> n_polyhedrons;
    m_polyhedrons.insert(m_polyhedrons.end(), n_polyhedrons, new Polyhedron);

    std::vector<std::size_t> polyhedrons_faces_nums;
    polyhedrons_faces_nums.reserve(n_polyhedrons);
    for (std::size_t i = 0; i < n_polyhedrons; i++) {
        std::size_t buf;
        input >> buf;
        polyhedrons_faces_nums.push_back(buf);
    }

    for (std::size_t i = 0; i < n_polyhedrons; i++) {
        for (std::size_t j = 0; j < polyhedrons_faces_nums[i]; j++) {
            std::size_t face_idx;
            input >> face_idx;

            for (std::size_t k = 0; k < 3; ++k) {
                if (!m_polyhedrons[i]->shell_contains(m_shell_faces[face_idx]->edges[k])) {
                    m_polyhedrons[i]->add_to_shell(m_shell_faces[face_idx]->edges[k]);

                    if (!m_polyhedrons[i]->shell_contains(m_shell_faces[face_idx]->edges[k]->verts[0]))
                        m_polyhedrons[i]->add_to_shell(m_shell_faces[face_idx]->edges[k]->verts[0]);

                    if (!m_polyhedrons[i]->shell_contains(m_shell_faces[face_idx]->edges[k]->verts[1]))
                        m_polyhedrons[i]->add_to_shell(m_shell_faces[face_idx]->edges[k]->verts[1]);
                }
            }

            m_polyhedrons[i]->add_to_shell(m_shell_faces[face_idx]);
        }
    }

    input.close();
}


void PolyhedralSet::set_poly_shell(const psg::PolyShell& poly_struct) {
    for (std::size_t i = 0; i < poly_struct.verts.size(); i++) {
        vec3 v(poly_struct.verts[i][0],
               poly_struct.verts[i][1],
               poly_struct.verts[i][2]);

        m_shell_verts.push_back(new shell::Vert(v));
    }

    for (const auto& face : poly_struct.faces) {
        if (!find_shell_edge(m_shell_verts[face[0]], m_shell_verts[face[1]])) {
            m_shell_edges.push_back(new shell::Edge(m_shell_verts[face[0]], m_shell_verts[face[1]]));
        }

        if (!find_shell_edge(m_shell_verts[face[1]], m_shell_verts[face[2]])) {
            m_shell_edges.push_back(new shell::Edge(m_shell_verts[face[1]], m_shell_verts[face[2]]));
        }

        if (!find_shell_edge(m_shell_verts[face[2]], m_shell_verts[face[0]])) {
            m_shell_edges.push_back(new shell::Edge(m_shell_verts[face[2]], m_shell_verts[face[0]]));
        }

        m_shell_faces.push_back(new shell::Face(
            find_shell_edge(m_shell_verts[face[0]], m_shell_verts[face[1]]),
            find_shell_edge(m_shell_verts[face[1]], m_shell_verts[face[2]]),
            find_shell_edge(m_shell_verts[face[2]], m_shell_verts[face[0]])));
    }

    m_polyhedrons.reserve(poly_struct.polyhs.size());
    for (std::size_t i = 0; i < poly_struct.polyhs.size(); ++i)
        m_polyhedrons.push_back(new Polyhedron(this));

    for (std::size_t i = 0; i < poly_struct.polyhs.size(); ++i) {
        for (std::size_t j = 0; j < poly_struct.polyhs[i].size(); ++j) {
            std::size_t face_idx = poly_struct.polyhs[i][j];
            for (std::size_t k = 0; k < 3; ++k) {
                if (!m_polyhedrons[i]->shell_contains(m_shell_faces[face_idx]->edges[k])) {
                    m_polyhedrons[i]->add_to_shell(m_shell_faces[face_idx]->edges[k]);

                    if (!m_polyhedrons[i]->shell_contains(m_shell_faces[face_idx]->edges[k]->verts[0]))
                        m_polyhedrons[i]->add_to_shell(m_shell_faces[face_idx]->edges[k]->verts[0]);

                    if (!m_polyhedrons[i]->shell_contains(m_shell_faces[face_idx]->edges[k]->verts[1]))
                        m_polyhedrons[i]->add_to_shell(m_shell_faces[face_idx]->edges[k]->verts[1]);
                }
            }

            m_polyhedrons[i]->add_to_shell(m_shell_faces[face_idx]);
        }
    }
}


std::string PolyhedralSet::output_filename(filetype filetype, std::string_view filename) const {
    if (filename != "_AUTO_")
        return filename.data();

    std::stringstream ss;
    ss << "phset_" << m_polyhedrons.size() << "_nph_";
    std::size_t av_nfe = 0;
    for (auto& polyhedr : m_polyhedrons)
        av_nfe += polyhedr->inner_tetrs().size();
    av_nfe /= m_polyhedrons.size();
    ss << av_nfe << "_phfe";
    switch (filetype) {
    case filetype::wavefront_obj:
        return ss.str() + ".obj";

    case filetype::lsdyna_keyword:
        return ss.str() + ".kw";
    }
    throw std::exception();
}


void PolyhedralSet::output(filetype filetype, std::string_view filename, unsigned polyhedral_set_id) {
    using namespace std::chrono;
    auto start = steady_clock::now();
    switch (filetype) {
    case filetype::wavefront_obj:
        output_obj(output_filename(filetype::wavefront_obj, filename));
        break;

    case filetype::lsdyna_keyword:
        output_lsdynakw(output_filename(filetype::lsdyna_keyword, filename), polyhedral_set_id);
        break;
    }
    auto stop = steady_clock::now();
    duration<double> elapsed = stop - start;

    m_log.mesh_file_writing_time = elapsed.count();
}




PolyhedralSet::PolyhedralSet(const psg::PolyShell& poly_struct) {
    set_poly_shell(poly_struct);
}


PolyhedralSet::~PolyhedralSet() {
    for (auto& polyhedr : m_polyhedrons)
        delete polyhedr;

    for (auto& face : m_shell_faces)
        delete face;
    for (auto& edge : m_shell_edges)
        delete edge;
    for (auto& vert : m_shell_verts)
        delete vert;
}




void PolyhedralSet::Log::write(std::ostream& stream) const {
    Logger logger(stream);
    logger << std::fixed;
    auto logWriteIfPositive = [&logger](std::string_view first, auto second, std::string_view third) {
        if constexpr (std::is_same<decltype(second), std::size_t>()) {
            logger << first.data() << second << third.data();
        } else {
            if (second >= 0.0)
                logger << first.data() << second << third.data();
            else logger << first.data() << "-" << "";
        }
    };

    logWriteIfPositive("Minimum quality", min_quality, "");
    logWriteIfPositive("Average quality", av_quality, "");
    logWriteIfPositive("Minimum mesh absGrad", min_mesh_abs_grad, "");
    logWriteIfPositive("Average mesh absGrad", av_mesh_abs_grad, "");
    logWriteIfPositive("Polyhedrons number", n_polyhs, "");
    logWriteIfPositive("Elements number", n_elems, "");
    logWriteIfPositive("Preferred edge length", pref_len, "");
    logWriteIfPositive("Shell triangulation time", shell_tr_time, "s");
    logWriteIfPositive("Volume exhaustion time", volume_exh_time, "s");
    logWriteIfPositive("Mesh file writing time", mesh_file_writing_time, "s");
}
