// Copyright © 2018-2019 Tokarev Artem Alekseevich. All rights reserved.
// Licensed under the MIT License.

#include "polyhedral-set.h"
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <ctime>
#include "helpers/spatial-algs/spatial-algs.h"

using namespace pmg;


//#define DEV_DEBUG




shell::Edge* PolyhedralSet::findShellEdge(const shell::Vert* v0, const shell::Vert* v1) const
{
    for (auto& s_edge : m_shellEdges)
    {
        if ((s_edge->verts[0] == v0  &&
             s_edge->verts[1] == v1) ||
            (s_edge->verts[1] == v0  &&
             s_edge->verts[0] == v1))
            return s_edge;
    }
    
    return nullptr;
}


void PolyhedralSet::triangulateShell()
{
    for (auto& svert : m_shellVerts)
        svert->attachedVert = new pmg::Vert(svert->pos());

    for (auto& sedge : m_shellEdges)
        sedge->segmentize(m_prefLen);

    size_t n_shell_faces = m_shellFaces.size();
    #pragma omp parallel for
    for (size_t i = 0; i < n_shell_faces; i++)
        m_shellFaces[i]->triangulate(m_prefLen);
}




void PolyhedralSet::outputObj(std::string_view filename) const
{
    std::ofstream file(filename.data());

    size_t index = 1;
    for (auto& svert : m_shellVerts)
    {
        file << "v " << (*svert)[0] << ' ' << (*svert)[1] << ' ' << (*svert)[2] << '\n';
        svert->attachedVert->globalNum = index++;
    }
    for (auto& sedge : m_shellEdges)
    {
        for (auto& vert : sedge->innerVerts())
        {
            file << "v " << (*vert)[0] << ' ' << (*vert)[1] << ' ' << (*vert)[2] << '\n';
            vert->globalNum = index++;
        }
    }
    for (auto& sface : m_shellFaces)
    {
        for (auto& vert : sface->innerVerts())
        {
            file << "v " << (*vert)[0] << ' ' << (*vert)[1] << ' ' << (*vert)[2] << '\n';
            vert->globalNum = index++;
        }
    }

    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& vert : polyhedr->innerVerts())
        {
            file << "v " << (*vert)[0] << ' ' << (*vert)[1] << ' ' << (*vert)[2] << '\n';
            vert->globalNum = index++;
        }
    }

    #ifdef DEV_DEBUG
    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& f_face : polyhedr->frontFaces())
        {
            auto face = f_face->face;

            std::vector<size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalNum) == gl_nums.end())
                        gl_nums.push_back(vert->globalNum);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
    #else
    for (auto& sface : m_shellFaces)
    {
        for (auto& face : sface->innerFaces())
        {
            std::vector<size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalNum) == gl_nums.end())
                        gl_nums.push_back(vert->globalNum);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }

    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& face : polyhedr->innerFaces())
        {
            std::vector<size_t> gl_nums;
            for (auto& edge : face->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalNum) == gl_nums.end())
                        gl_nums.push_back(vert->globalNum);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
    #endif

    file.close();
}


void PolyhedralSet::outputLSDynaKeyword_PART(std::ofstream& file) const
{
    for (size_t id = 1; id < m_polyhedrons.size() + 1; id++)
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
    size_t index = 1;
    for (auto& svert : m_shellVerts)
    {
        svert->attachedVert->globalNum = index++;
        file << std::setw(8)  << svert->attachedVert->globalNum;
        file << std::setw(16) << (*svert)[0];
        file << std::setw(16) << (*svert)[1];
        file << std::setw(16) << (*svert)[2];
        file << "       0       0\n";
    }
    for (auto& sedge : m_shellEdges)
    {
        for (auto& vert : sedge->innerVerts())
        {
            vert->globalNum = index++;
            file << std::setw(8)  << vert->globalNum;
            file << std::setw(16) << (*vert)[0];
            file << std::setw(16) << (*vert)[1];
            file << std::setw(16) << (*vert)[2];
            file << "       0       0\n";
        }
    }
    for (auto& sface : m_shellFaces)
    {
        for (auto& vert : sface->innerVerts())
        {
            vert->globalNum = index++;
            file << std::setw(8)  << vert->globalNum;
            file << std::setw(16) << (*vert)[0];
            file << std::setw(16) << (*vert)[1];
            file << std::setw(16) << (*vert)[2];
            file << "       0       0\n";
        }
    }

    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& vert : polyhedr->innerVerts())
        {
            vert->globalNum = index++;
            file << std::setw(8)  << vert->globalNum;
            file << std::setw(16) << (*vert)[0];
            file << std::setw(16) << (*vert)[1];
            file << std::setw(16) << (*vert)[2];
            file << "       0       0\n";
        }
    }
}


void PolyhedralSet::outputLSDynaKeyword_ELEMENT_SOLID(std::ofstream& file, unsigned PolyhedralSetId) const
{
    file << "*ELEMENT_SOLID\n"
        "$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n";
    size_t pid = 1;
    size_t eid = 10000000 * PolyhedralSetId + 1;
    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& tetr : polyhedr->innerTetrs())
        {
            file << std::setw(8) << eid++;
            file << std::setw(8) << pid;
            file << std::setw(8) << tetr->verts[0]->globalNum;
            vec3 v0 = tetr->verts[1]->pos() - tetr->verts[0]->pos();
            vec3 v1 = tetr->verts[2]->pos() - tetr->verts[0]->pos();
            vec3 v2 = tetr->verts[3]->pos() - tetr->verts[0]->pos();
            if (vec3::dot(v2, vec3::cross(v0, v1)) > static_cast<real_t>(0.0))
            {
                file << std::setw(8) << tetr->verts[1]->globalNum;
                file << std::setw(8) << tetr->verts[2]->globalNum;
            }
            else
            {
                file << std::setw(8) << tetr->verts[2]->globalNum;
                file << std::setw(8) << tetr->verts[1]->globalNum;
            }
            file << std::setw(8) << tetr->verts[3]->globalNum;
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
    size_t av_nfe = 0;
    for (auto& polyhedr : m_polyhedrons)
        av_nfe += polyhedr->innerTetrs().size();
    av_nfe /= m_polyhedrons.size();
    ss << av_nfe << "_phfe.log";

    return ss.str();
}


void PolyhedralSet::generateMesh(real_t preferredLength)
{
    m_prefLen = preferredLength;

    std::clock_t shell_triang_start = std::clock();
    try
    {
        triangulateShell();
    }
    catch (std::logic_error error)
    {
        output(FileType::WavefrontObj, "debug.obj");
        throw error;
    }
    double shell_triang_elapsed = static_cast<double>(std::clock() - shell_triang_start) / CLOCKS_PER_SEC;

    std::clock_t volume_exh_start = std::clock();
    size_t n_polyhs = m_polyhedrons.size();
    #pragma omp parallel for
    for (size_t i = 0; i < n_polyhs; i++)
    {
        try
        {
            m_polyhedrons[i]->generateMesh(m_prefLen);
        }
        catch (std::logic_error error)
        {
            output(FileType::WavefrontObj, "debug.obj");
            throw error;
        }
    }
    double volume_exh_elapsed = static_cast<double>(std::clock() - volume_exh_start) / CLOCKS_PER_SEC;

    m_log.nPolyhs = m_polyhedrons.size();
    m_log.prefLen = m_prefLen;
    m_log.shellTrTime   = shell_triang_elapsed;
    m_log.volumeExhTime = volume_exh_elapsed;
    m_isLogged = false;
}


void PolyhedralSet::optimizeMesh(settings::Optimization optSettings)
{
    size_t n_polyhs = m_polyhedrons.size();
    #pragma omp parallel for
    for (size_t i = 0; i < n_polyhs; i++)
        m_polyhedrons[i]->optimizeMesh(optSettings);

    m_isLogged = false;
}


const PolyMesh* PolyhedralSet::structurizeMesh()
{
    m_lastMesh = new PolyMesh;

    m_lastMesh->nCryses = m_polyhedrons.size();
    m_lastMesh->nCrysesTetrs = new size_t[m_lastMesh->nCryses];

    size_t nodes_num = m_shellVerts.size();
    for (auto& sedge : m_shellEdges)
        nodes_num += sedge->innerVerts().size();
    for (auto& sface : m_shellFaces)
        nodes_num += sface->innerVerts().size();

    for (auto& polyhedr : m_polyhedrons)
        nodes_num += static_cast<size_t>(std::count_if(
            polyhedr->innerVerts().begin(),
            polyhedr->innerVerts().end(),
            [](auto vert) { return static_cast<bool>(vert); }));

    m_lastMesh->nNodes = nodes_num;
    m_lastMesh->nodesPositions = new real_t[3 * m_lastMesh->nNodes];
    Vert** verts_ptrs = new Vert*[m_lastMesh->nNodes];

    auto lastMesh_buf = m_lastMesh;
    auto FindNodeIndex = [lastMesh_buf, verts_ptrs](Vert* node_ptr)-> size_t
    {
        for (size_t i = 0; i < lastMesh_buf->nNodes; i++)
            if (verts_ptrs[i] == node_ptr)
                return i;

        return lastMesh_buf->nNodes;
    };

    size_t vert_ind = 0;
    for (auto& svert : m_shellVerts)
    {
        verts_ptrs[vert_ind] = svert->attachedVert;
        m_lastMesh->nodesPositions[3 * vert_ind]     = (*svert)[0];
        m_lastMesh->nodesPositions[3 * vert_ind + 1] = (*svert)[1];
        m_lastMesh->nodesPositions[3 * vert_ind + 2] = (*svert)[2];
        vert_ind++;
    }
    for (auto& sedge : m_shellEdges)
    {
        for (auto& vert : sedge->innerVerts())
        {
            verts_ptrs[vert_ind] = vert;
            m_lastMesh->nodesPositions[3 * vert_ind]     = (*vert)[0];
            m_lastMesh->nodesPositions[3 * vert_ind + 1] = (*vert)[1];
            m_lastMesh->nodesPositions[3 * vert_ind + 2] = (*vert)[2];
            vert_ind++;
        }
    }
    for (auto& sface : m_shellFaces)
    {
        for (auto& vert : sface->innerVerts())
        {
            verts_ptrs[vert_ind] = vert;
            m_lastMesh->nodesPositions[3 * vert_ind]     = (*vert)[0];
            m_lastMesh->nodesPositions[3 * vert_ind + 1] = (*vert)[1];
            m_lastMesh->nodesPositions[3 * vert_ind + 2] = (*vert)[2];
            vert_ind++;
        }
    }
    for (auto& polyhedr : m_polyhedrons)
    {
        for (auto& vert : polyhedr->innerVerts())
        {
            verts_ptrs[vert_ind] = vert;
            m_lastMesh->nodesPositions[3 * vert_ind]     = (*vert)[0];
            m_lastMesh->nodesPositions[3 * vert_ind + 1] = (*vert)[1];
            m_lastMesh->nodesPositions[3 * vert_ind + 2] = (*vert)[2];
            vert_ind++;
        }
    }

    size_t tetrs_num = 0;
    for (size_t i = 0; i < m_lastMesh->nCryses; i++)
    {
        m_lastMesh->nCrysesTetrs[i] = static_cast<size_t>(std::count_if(
            m_polyhedrons[i]->innerTetrs().begin(),
            m_polyhedrons[i]->innerTetrs().end(),
            [](auto tetr) { return static_cast<bool>(tetr); }));

        tetrs_num += m_lastMesh->nCrysesTetrs[i];
    }

    m_lastMesh->nTetrs = tetrs_num;
    m_lastMesh->tetrs = new size_t[4 * m_lastMesh->nTetrs];

    size_t tetr_ind = 0;
    for (size_t i = 0; i < m_lastMesh->nCryses; i++)
    {
        for (auto& tetr : m_polyhedrons[i]->innerTetrs())
        {
            m_lastMesh->tetrs[4 * tetr_ind]     = FindNodeIndex(tetr->verts[0]);
            m_lastMesh->tetrs[4 * tetr_ind + 1] = FindNodeIndex(tetr->verts[1]);
            m_lastMesh->tetrs[4 * tetr_ind + 2] = FindNodeIndex(tetr->verts[1]);
            m_lastMesh->tetrs[4 * tetr_ind + 3] = FindNodeIndex(tetr->verts[1]);
            tetr_ind++;
        }
    }

    return m_lastMesh;
}


const PolyMesh* PolyhedralSet::getLastMesh()
{
    return m_lastMesh;
}




PolyhedralSet::Log PolyhedralSet::log()
{
    if (m_isLogged)
        return m_log;

    real_t min_q = 1.0, av_q = 0.0;
    size_t n_elems = 0;
    for (auto& polyh : m_polyhedrons)
    {
        auto buf_min_av_q = polyh->analyzeMeshQuality();
        if (buf_min_av_q.first < min_q)
            min_q = buf_min_av_q.first;
        av_q += buf_min_av_q.second;
        n_elems += polyh->innerTetrs().size();
    }
    av_q /= m_polyhedrons.size();

    m_log.minMeshAbsGrad;
    m_log.avMeshAbsGrad;
    m_log.minQuality = min_q;
    m_log.avQuality  = av_q;
    m_log.nElems     = n_elems;

    return m_log;
}




void PolyhedralSet::input(std::string_view polyStructFileName)
{
    std::ifstream input(polyStructFileName.data());

    size_t nodes_num;
    input >> nodes_num;

    for (size_t i = 0; i < nodes_num; i++)
    {
        real_t coors[3];
        input >> coors[0];
        input >> coors[1];
        input >> coors[2];

        m_shellVerts.push_back(new shell::Vert(coors[0], coors[1], coors[2]));
    }

    size_t faces_num;
    input >> faces_num;

    for (size_t i = 0; i < faces_num; i++)
    {
        size_t face_nodes_inds[3];
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

    size_t polyhedrons_num;
    input >> polyhedrons_num;
    m_polyhedrons.insert(m_polyhedrons.end(), polyhedrons_num, new Polyhedron);

    size_t* polyhedrons_faces_nums = new size_t[polyhedrons_num];
    for (size_t i = 0; i < polyhedrons_num; i++)
        input >> polyhedrons_faces_nums[i];

    for (size_t i = 0; i < polyhedrons_num; i++)
    {
        for (size_t j = 0; j < polyhedrons_faces_nums[i]; j++)
        {
            size_t face_ind;
            input >> face_ind;

            for (int k = 0; k < 3; ++k)
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


void PolyhedralSet::input(const psg::PolyShell& polyStruct)
{
    for (size_t i = 0; i < polyStruct.verts.size(); i++)
    {
        real_t coors[3];
        coors[0] = polyStruct.verts[i][0];
        coors[1] = polyStruct.verts[i][1];
        coors[2] = polyStruct.verts[i][2];

        m_shellVerts.push_back(new shell::Vert(coors[0], coors[1], coors[2]));
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

    m_polyhedrons.reserve(polyStruct.polyhedrons.size());
    for (size_t i = 0; i < polyStruct.polyhedrons.size(); ++i)
        m_polyhedrons.push_back(new Polyhedron(this));

    for (size_t i = 0; i < polyStruct.polyhedrons.size(); ++i)
    {
        for (size_t j = 0; j < polyStruct.polyhedrons[i].size(); ++j)
        {
            size_t face_ind = polyStruct.polyhedrons[i][j];

            for (int k = 0; k < 3; ++k)
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
}


std::string PolyhedralSet::generateOutputFilename(FileType filetype, std::string_view filename) const
{
    if (filename != "_AUTO_")
        return filename.data();

    std::stringstream ss;
    ss << "phset_" << m_polyhedrons.size() << "_nph_";
    size_t av_nfe = 0;
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


void PolyhedralSet::output(FileType filetype, std::string_view filename, unsigned PolyhedralSetId)
{
    std::clock_t start = std::clock();
    switch (filetype)
    {
    case FileType::WavefrontObj:
        outputObj(generateOutputFilename(FileType::WavefrontObj, filename));
        break;

    case FileType::LsDynaKeyword:
        outputLSDynaKeyword(generateOutputFilename(FileType::LsDynaKeyword, filename), PolyhedralSetId);
        break;
    }
    double elapsed = static_cast<double>(std::clock() - start) / CLOCKS_PER_SEC;

    m_log.meshFileWritingTime = elapsed;
}




PolyhedralSet::PolyhedralSet() {}


PolyhedralSet::PolyhedralSet(std::string_view polyStructFileName)
{
    input(polyStructFileName);
}


PolyhedralSet::PolyhedralSet(const psg::PolyShell& polyStruct)
{
    input(polyStruct);
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
        if constexpr (std::is_same<decltype(second), size_t>())
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
    logWriteIfPositive("Minimum absMeshGrad",      minMeshAbsGrad, "");
    logWriteIfPositive("Average absMeshGrad",      avMeshAbsGrad,  "");
    logWriteIfPositive("Polyhedrons number",       nPolyhs,        "");
    logWriteIfPositive("Elements number",          nElems,         "");
    logWriteIfPositive("Preferred edge length",    prefLen,        "");
    logWriteIfPositive("Shell triangulation time", shellTrTime,         "s");
    logWriteIfPositive("Volume exhaustion time",   volumeExhTime,       "s");
    logWriteIfPositive("Optimization time",        optimizationTime,    "s");
    logWriteIfPositive("Mesh file writing time",   meshFileWritingTime, "s");
}
