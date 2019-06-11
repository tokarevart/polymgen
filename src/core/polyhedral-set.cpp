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
using spt::vec3;


//#define DEV_DEBUG


shell::Edge* PolyhedralSet::findShellEdge(const shell::Vert* v0, const shell::Vert* v1) const
{
    for (auto& s_edge : m_shellEdges)
        if (   (   s_edge->verts[0] == v0
                && s_edge->verts[1] == v1)
            || (   s_edge->verts[1] == v0
                && s_edge->verts[0] == v1))
            return s_edge;
    
    return nullptr;
}


void PolyhedralSet::triangulateShell(genparams::Shell genParams)
{
    for (auto& svert : m_shellVerts)
        svert->attachedVert = new pmg::Vert(svert->pos());

    for (auto& sedge : m_shellEdges)
        sedge->segmentize(m_prefLen);

    std::size_t n_shell_faces = m_shellFaces.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_shell_faces; i++)
        m_shellFaces[i]->triangulate(m_prefLen, genParams);
}


void PolyhedralSet::outputObj(std::string_view filename) const
{
    std::ofstream file(filename.data());

    std::size_t index = 1;
    for (auto& svert : m_shellVerts)
    {
        file << "v " << svert->pos()[0] << ' ' << svert->pos()[1] << ' ' << svert->pos()[2] << '\n';
        svert->attachedVert->globalIdx = index++;
    }
    for (auto& sedge : m_shellEdges)
    {
        for (auto& vert : sedge->innerVerts())
        {
            file << "v " << vert->pos()[0] << ' ' << vert->pos()[1] << ' ' << vert->pos()[2] << '\n';
            vert->globalIdx = index++;
        }
    }
    for (auto& sface : m_shellFaces)
    {
        for (auto& vert : sface->innerVerts())
        {
            file << "v " << vert->pos()[0] << ' ' << vert->pos()[1] << ' ' << vert->pos()[2] << '\n';
            vert->globalIdx = index++;
        }
    }

    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& vert : polyhedr->innerVerts())
        {
            file << "v " << vert->pos()[0] << ' ' << vert->pos()[1] << ' ' << vert->pos()[2] << '\n';
            vert->globalIdx = index++;
        }
    }

    #ifdef DEV_DEBUG
    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& f_face : polyhedr->frontFaces())
        {
            auto face = f_face->face;

            std::vector<std::size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalIdx) == gl_nums.end())
                        gl_nums.push_back(vert->globalIdx);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
    #else
    for (auto& sface : m_shellFaces)
    {
        for (auto& face : sface->innerFaces())
        {
            std::vector<std::size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalIdx) == gl_nums.end())
                        gl_nums.push_back(vert->globalIdx);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }

    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& face : polyhedr->innerFaces())
        {
            std::vector<std::size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalIdx) == gl_nums.end())
                        gl_nums.push_back(vert->globalIdx);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
    #endif

    file.close();
}


void PolyhedralSet::outputLSDynaKeyword_PART(std::ofstream& file) const
{
    for (std::size_t id = 1; id < m_polyhedrons.size() + 1; id++)
    {
        file << "*PART\n"
            "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n";
        file << std::setw(10) << id;
        file << "         0";
        file << std::setw(10) << id;
        file << "         0         0         0         0         0\n";
    }
}


void PolyhedralSet::outputLSDynaKeyword_NODE(std::ofstream& file) const
{
    file << "*NODE\n"
        "$#   nid               x               y               z      tc      rc\n";
    std::size_t index = 1;
    for (auto& svert : m_shellVerts)
    {
        svert->attachedVert->globalIdx = index++;
        file << std::setw(8)  << svert->attachedVert->globalIdx;
        file << std::setw(16) << svert->pos()[0];
        file << std::setw(16) << svert->pos()[1];
        file << std::setw(16) << svert->pos()[2];
        file << "       0       0\n";
    }
    for (auto& sedge : m_shellEdges)
    {
        for (auto& vert : sedge->innerVerts())
        {
            vert->globalIdx = index++;
            file << std::setw(8)  << vert->globalIdx;
            file << std::setw(16) << vert->pos()[0];
            file << std::setw(16) << vert->pos()[1];
            file << std::setw(16) << vert->pos()[2];
            file << "       0       0\n";
        }
    }
    for (auto& sface : m_shellFaces)
    {
        for (auto& vert : sface->innerVerts())
        {
            vert->globalIdx = index++;
            file << std::setw(8)  << vert->globalIdx;
            file << std::setw(16) << vert->pos()[0];
            file << std::setw(16) << vert->pos()[1];
            file << std::setw(16) << vert->pos()[2];
            file << "       0       0\n";
        }
    }

    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& vert : polyhedr->innerVerts())
        {
            vert->globalIdx = index++;
            file << std::setw(8)  << vert->globalIdx;
            file << std::setw(16) << vert->pos()[0];
            file << std::setw(16) << vert->pos()[1];
            file << std::setw(16) << vert->pos()[2];
            file << "       0       0\n";
        }
    }
}


void PolyhedralSet::outputLSDynaKeyword_ELEMENT_SOLID(std::ofstream& file, unsigned PolyhedralSetId) const
{
    file << "*ELEMENT_SOLID\n"
        "$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n";
    std::size_t pid = 1;
    std::size_t eid = 10000000 * PolyhedralSetId + 1;
    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& tetr : polyhedr->innerTetrs())
        {
            file << std::setw(8) << eid++;
            file << std::setw(8) << pid;
            file << std::setw(8) << tetr->verts[0]->globalIdx;
            vec3 v0 = tetr->verts[1]->pos() - tetr->verts[0]->pos();
            vec3 v1 = tetr->verts[2]->pos() - tetr->verts[0]->pos();
            vec3 v2 = tetr->verts[3]->pos() - tetr->verts[0]->pos();
            if (vec3::dot(v2, vec3::cross(v0, v1)) > static_cast<real_t>(0.0))
            {
                file << std::setw(8) << tetr->verts[1]->globalIdx;
                file << std::setw(8) << tetr->verts[2]->globalIdx;
            }
            else
            {
                file << std::setw(8) << tetr->verts[2]->globalIdx;
                file << std::setw(8) << tetr->verts[1]->globalIdx;
            }
            file << std::setw(8) << tetr->verts[3]->globalIdx;
            file << "       0       0       0       0\n";
        }
        pid++;
    }
}


void PolyhedralSet::outputLSDynaKeyword(const std::string& filename, unsigned PolyhedralSetId) const
{
    std::ofstream file(filename);
    file.setf(std::ios::right | std::ios::fixed);

    outputLSDynaKeyword_PART(file);
    outputLSDynaKeyword_NODE(file);
    outputLSDynaKeyword_ELEMENT_SOLID(file, PolyhedralSetId);
    file << "*END";

    file.close();
}


std::string PolyhedralSet::generateLogFileName() const
{
    std::stringstream ss;
    ss << "pmg_log_" << m_polyhedrons.size() << "_nph_";
    std::size_t av_nfe = 0;
    for (auto& polyhedr : m_polyhedrons)
        av_nfe += polyhedr->innerTetrs().size();
    av_nfe /= m_polyhedrons.size();
    ss << av_nfe << "_phfe.log";

    return ss.str();
}


void PolyhedralSet::tetrahedralize(real_t preferredLength, genparams::Polyhedron genParams)
{
    m_prefLen = preferredLength;

    auto shell_triang_start = std::chrono::steady_clock::now();
    try
    {
        triangulateShell(genParams.shell);
    }
    catch (std::logic_error error)
    {
        output(FileType::WavefrontObj, "debug.obj");
        throw error;
    }
    auto shell_triang_stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> shell_triang_elapsed = shell_triang_stop - shell_triang_start;

    auto volume_exh_start = std::chrono::steady_clock::now();
    std::size_t n_polyhs = m_polyhedrons.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_polyhs; i++)
    {
        try
        {
            m_polyhedrons[i]->tetrahedralize(m_prefLen, genParams.volume);
        }
        catch (std::logic_error error)
        {
            output(FileType::WavefrontObj, "debug.obj");
            throw error;
        }
    }
    auto volume_exh_stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> volume_exh_elapsed = volume_exh_stop - volume_exh_start;

    m_log.nPolyhs = m_polyhedrons.size();
    m_log.prefLen = m_prefLen;
    m_log.shellTrTime   = shell_triang_elapsed.count();
    m_log.volumeExhTime = volume_exh_elapsed.count();
    m_isLogged = false;
}


void PolyhedralSet::smoothMesh(std::size_t nItersVolume, std::size_t nItersShell)
{
    std::size_t n_sfaces = m_shellFaces.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_sfaces; i++)
        m_shellFaces[i]->smoothMesh(nItersShell);

    std::size_t n_polyhs = m_polyhedrons.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_polyhs; i++)
        m_polyhedrons[i]->smoothMesh(nItersVolume);

    m_isLogged = false;
}


void PolyhedralSet::shellDelaunayPostP()
{
    std::size_t n_sfaces = m_shellFaces.size();
    #pragma omp parallel for
    for (std::size_t i = 0; i < n_sfaces; i++)
        m_shellFaces[i]->delaunayPostP();

    m_isLogged = false;
}


PolyMesh PolyhedralSet::mesh() const
{
    PolyMesh mesh;
    mesh.polyhs.assign(m_polyhedrons.size(), std::vector<PolyMesh::TetrIdx>());

    std::size_t n_verts = m_shellVerts.size();
    for (auto& sedge : m_shellEdges)
        n_verts += sedge->innerVerts().size();
    for (auto& sface : m_shellFaces)
        n_verts += sface->innerVerts().size();
    for (auto& polyh : m_polyhedrons)
        n_verts += static_cast<std::size_t>(std::count_if(
            polyh->innerVerts().begin(),
            polyh->innerVerts().end(),
            [](auto vert) { return static_cast<bool>(vert); }));

    mesh.verts.assign(n_verts, PolyMesh::Vert());

    std::size_t vert_idx = 0;
    for (auto& svert : m_shellVerts)
    {
        mesh.verts[vert_idx] = svert->pos().x;
        svert->attachedVert->globalIdx = vert_idx;
        vert_idx++;
    }
    auto addVertsToLastMesh = [&mesh, &vert_idx](auto verts) mutable
    {
        for (auto& vert : verts)
        {
            mesh.verts[vert_idx] = vert->pos().x;
            vert->globalIdx = vert_idx;
            vert_idx++;
        }
    };
    for (auto& sedge : m_shellEdges)  addVertsToLastMesh(sedge->innerVerts());
    for (auto& sface : m_shellFaces)  addVertsToLastMesh(sface->innerVerts());
    for (auto& polyh : m_polyhedrons) addVertsToLastMesh(polyh->innerVerts());

    std::size_t n_tetrs = 0;
    for (std::size_t i = 0; i < mesh.polyhs.size(); i++)
    {
        mesh.polyhs[i].assign(m_polyhedrons[i]->innerTetrs().size(), 0);
        n_tetrs += mesh.polyhs[i].size();
    }

    mesh.tetrs.assign(n_tetrs, { 0, 0, 0, 0 });

    std::size_t tetr_idx = 0;
    for (auto& polyh : m_polyhedrons)
        for (auto& tetr : polyh->innerTetrs())
        {
            std::array<PolyMesh::VertIdx, 4> verts_idces;
            verts_idces = { tetr->verts[0]->globalIdx,
                            tetr->verts[1]->globalIdx,
                            tetr->verts[2]->globalIdx,
                            tetr->verts[3]->globalIdx };

            vec3 v0 = tetr->verts[1]->pos() - tetr->verts[0]->pos();
            vec3 v1 = tetr->verts[2]->pos() - tetr->verts[0]->pos();
            vec3 v2 = tetr->verts[3]->pos() - tetr->verts[0]->pos();
            if (vec3::dot(v2, vec3::cross(v0, v1)) < static_cast<real_t>(0.0))
                std::swap(verts_idces[0], verts_idces[1]);

            mesh.tetrs[tetr_idx++] = verts_idces;
        }

    return mesh;
}


PolyhedralSet::Log PolyhedralSet::log()
{
    if (m_isLogged)
        return m_log;

    std::size_t n_elems = 0;
    real_t min_q = 1.0, av_q = 0.0;
    real_t min_g = 1.0, av_g = 0.0;
    for (auto& polyh : m_polyhedrons)
    {
        auto cur_min_av_q = polyh->analyzeMeshQuality();
        auto cur_min_av_g = polyh->analyzeMeshAbsGrad();
        min_q = std::min(cur_min_av_q.first, min_q);
        min_g = std::min(cur_min_av_g.first, min_g);
        av_q += cur_min_av_q.second;
        av_g += cur_min_av_g.second;
        n_elems += polyh->innerTetrs().size();
    }
    av_q /= m_polyhedrons.size();
    av_g /= m_polyhedrons.size();

    m_log.minQuality = min_q;
    m_log.avQuality  = av_q;
    m_log.minMeshAbsGrad = min_g;
    m_log.avMeshAbsGrad  = av_g;
    m_log.nElems = n_elems;

    return m_log;
}


void PolyhedralSet::setPolyShell(std::string_view polyStructFileName)
{
    std::ifstream input(polyStructFileName.data());

    std::size_t nodes_num;
    input >> nodes_num;

    for (std::size_t i = 0; i < nodes_num; i++)
    {
        spt::vec3 v;
        input >> v.x[0];
        input >> v.x[1];
        input >> v.x[2];

        m_shellVerts.push_back(new shell::Vert(v));
    }

    std::size_t faces_num;
    input >> faces_num;

    for (std::size_t i = 0; i < faces_num; i++)
    {
        std::array<std::size_t, 3> face_nodes_inds;
        input >> face_nodes_inds[0];
        input >> face_nodes_inds[1];
        input >> face_nodes_inds[2];

        if (!findShellEdge(m_shellVerts[face_nodes_inds[0]], m_shellVerts[face_nodes_inds[1]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[face_nodes_inds[0]], m_shellVerts[face_nodes_inds[1]]));
        }

        if (!findShellEdge(m_shellVerts[face_nodes_inds[1]], m_shellVerts[face_nodes_inds[2]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[face_nodes_inds[1]], m_shellVerts[face_nodes_inds[2]]));
        }

        if (!findShellEdge(m_shellVerts[face_nodes_inds[2]], m_shellVerts[face_nodes_inds[0]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[face_nodes_inds[2]], m_shellVerts[face_nodes_inds[0]]));
        }

        m_shellFaces.push_back(new shell::Face(
            findShellEdge(m_shellVerts[face_nodes_inds[0]], m_shellVerts[face_nodes_inds[1]]),
            findShellEdge(m_shellVerts[face_nodes_inds[1]], m_shellVerts[face_nodes_inds[2]]),
            findShellEdge(m_shellVerts[face_nodes_inds[2]], m_shellVerts[face_nodes_inds[0]])));
    }

    std::size_t polyhedrons_num;
    input >> polyhedrons_num;
    m_polyhedrons.insert(m_polyhedrons.end(), polyhedrons_num, new Polyhedron);

    std::size_t* polyhedrons_faces_nums = new std::size_t[polyhedrons_num];
    for (std::size_t i = 0; i < polyhedrons_num; i++)
        input >> polyhedrons_faces_nums[i];

    for (std::size_t i = 0; i < polyhedrons_num; i++)
    {
        for (std::size_t j = 0; j < polyhedrons_faces_nums[i]; j++)
        {
            std::size_t face_ind;
            input >> face_ind;

            for (std::size_t k = 0; k < 3; ++k)
            {
                if (!m_polyhedrons[i]->shellContains(m_shellFaces[face_ind]->edges[k]))
                {
                    m_polyhedrons[i]->addToShell(m_shellFaces[face_ind]->edges[k]);

                    if (!m_polyhedrons[i]->shellContains(m_shellFaces[face_ind]->edges[k]->verts[0]))
                        m_polyhedrons[i]->addToShell(m_shellFaces[face_ind]->edges[k]->verts[0]);

                    if (!m_polyhedrons[i]->shellContains(m_shellFaces[face_ind]->edges[k]->verts[1]))
                        m_polyhedrons[i]->addToShell(m_shellFaces[face_ind]->edges[k]->verts[1]);
                }
            }

            m_polyhedrons[i]->addToShell(m_shellFaces[face_ind]);
        }
    }

    input.close();
    delete[] polyhedrons_faces_nums;
}


void PolyhedralSet::setPolyShell(const psg::PolyShell& polyStruct)
{
    for (std::size_t i = 0; i < polyStruct.verts.size(); i++)
    {
        spt::vec3 v(polyStruct.verts[i][0],
                    polyStruct.verts[i][1],
                    polyStruct.verts[i][2]);

        m_shellVerts.push_back(new shell::Vert(v));
    }

    for (const auto& face : polyStruct.faces)
    {
        if (!findShellEdge(m_shellVerts[face[0]], m_shellVerts[face[1]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[face[0]], m_shellVerts[face[1]]));
        }

        if (!findShellEdge(m_shellVerts[face[1]], m_shellVerts[face[2]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[face[1]], m_shellVerts[face[2]]));
        }

        if (!findShellEdge(m_shellVerts[face[2]], m_shellVerts[face[0]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[face[2]], m_shellVerts[face[0]]));
        }

        m_shellFaces.push_back(new shell::Face(
            findShellEdge(m_shellVerts[face[0]], m_shellVerts[face[1]]),
            findShellEdge(m_shellVerts[face[1]], m_shellVerts[face[2]]),
            findShellEdge(m_shellVerts[face[2]], m_shellVerts[face[0]])));
    }

    m_polyhedrons.reserve(polyStruct.polyhs.size());
    for (std::size_t i = 0; i < polyStruct.polyhs.size(); ++i)
        m_polyhedrons.push_back(new Polyhedron(this));

    for (std::size_t i = 0; i < polyStruct.polyhs.size(); ++i)
    {
        for (std::size_t j = 0; j < polyStruct.polyhs[i].size(); ++j)
        {
            std::size_t face_idx = polyStruct.polyhs[i][j];
            for (std::size_t k = 0; k < 3; ++k)
            {
                if (!m_polyhedrons[i]->shellContains(m_shellFaces[face_idx]->edges[k]))
                {
                    m_polyhedrons[i]->addToShell(m_shellFaces[face_idx]->edges[k]);

                    if (!m_polyhedrons[i]->shellContains(m_shellFaces[face_idx]->edges[k]->verts[0]))
                        m_polyhedrons[i]->addToShell(m_shellFaces[face_idx]->edges[k]->verts[0]);

                    if (!m_polyhedrons[i]->shellContains(m_shellFaces[face_idx]->edges[k]->verts[1]))
                        m_polyhedrons[i]->addToShell(m_shellFaces[face_idx]->edges[k]->verts[1]);
                }
            }

            m_polyhedrons[i]->addToShell(m_shellFaces[face_idx]);
        }
    }
}


std::string PolyhedralSet::generateOutputFilename(FileType filetype, std::string_view filename) const
{
    if (filename != "_AUTO_")
        return filename.data();

    std::stringstream ss;
    ss << "phset_" << m_polyhedrons.size() << "_nph_";
    std::size_t av_nfe = 0;
    for (auto& polyhedr : m_polyhedrons)
        av_nfe += polyhedr->innerTetrs().size();
    av_nfe /= m_polyhedrons.size();
    ss << av_nfe << "_phfe";
    switch (filetype)
    {
    case FileType::WavefrontObj:
        return ss.str() + ".obj";

    case FileType::LsDynaKeyword:
        return ss.str() + ".kw";
    }
    throw std::exception();
}


void PolyhedralSet::output(FileType filetype, std::string_view filename, unsigned polyhedralSetId)
{
    auto start = std::chrono::steady_clock::now();
    switch (filetype)
    {
    case FileType::WavefrontObj:
        outputObj(generateOutputFilename(FileType::WavefrontObj, filename));
        break;

    case FileType::LsDynaKeyword:
        outputLSDynaKeyword(generateOutputFilename(FileType::LsDynaKeyword, filename), polyhedralSetId);
        break;
    }
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;

    m_log.meshFileWritingTime = elapsed.count();
}




PolyhedralSet::PolyhedralSet(const psg::PolyShell& polyStruct)
{
    setPolyShell(polyStruct);
}


PolyhedralSet::~PolyhedralSet()
{
    for (auto& polyhedr : m_polyhedrons)
        delete polyhedr;

    for (auto& face : m_shellFaces)
        delete face;
    for (auto& edge : m_shellEdges)
        delete edge;
    for (auto& vert : m_shellVerts)
        delete vert;
}




void PolyhedralSet::Log::write(std::ostream& stream) const
{
    Logger logger(stream);
    logger << std::fixed;
    auto logWriteIfPositive = [&logger](std::string_view first, auto second, std::string_view third)
    {
        if constexpr (std::is_same<decltype(second), std::size_t>())
        {
             logger << first.data() << second << third.data();
        }
        else
        {
            if (second >= 0.0)
                 logger << first.data() << second << third.data();
            else logger << first.data() << "-" << "";
        }
    };

    logWriteIfPositive("Minimum quality",          minQuality,     "");
    logWriteIfPositive("Average quality",          avQuality,      "");
    logWriteIfPositive("Minimum mesh absGrad",     minMeshAbsGrad, "");
    logWriteIfPositive("Average mesh absGrad",     avMeshAbsGrad,  "");
    logWriteIfPositive("Polyhedrons number",       nPolyhs,        "");
    logWriteIfPositive("Elements number",          nElems,         "");
    logWriteIfPositive("Preferred edge length",    prefLen,        "");
    logWriteIfPositive("Shell triangulation time", shellTrTime,         "s");
    logWriteIfPositive("Volume exhaustion time",   volumeExhTime,       "s");
    logWriteIfPositive("Mesh file writing time",   meshFileWritingTime, "s");
}
