#include "polycrystal.h"
#include <stddef.h>
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "helpers/spatial-algs/spatial-algs.h"
#include "helpers/timer.h"

using namespace pmg;


//#define DEV_DEBUG




shell::Edge* Polycrystal::findShellEdge(const shell::Vertex* v0, const shell::Vertex* v1) const
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


void Polycrystal::triangulateShell()
{
    for (auto& svert : m_shellVerts)
        svert->attachedVert = new pmg::Vertex(svert->pos());

    for (auto& sedge : m_shellEdges)
        sedge->segmentize(m_preferredLength);

    size_t n_shell_facets = m_shellFacets.size();
    #pragma omp parallel for
    for (size_t i = 0; i < n_shell_facets; i++)
        m_shellFacets[i]->triangulate(m_preferredLength);
}




void Polycrystal::outputObj(const std::string& filename) const
{
    std::ofstream file(filename);

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
    for (auto& sfacet : m_shellFacets)
    {
        for (auto& vert : sfacet->innerVerts())
        {
            file << "v " << (*vert)[0] << ' ' << (*vert)[1] << ' ' << (*vert)[2] << '\n';
            vert->globalNum = index++;
        }
    }

    for (auto& crys : m_crystallites)
    {
        for (auto& vert : crys->innerVerts())
        {
            file << "v " << (*vert)[0] << ' ' << (*vert)[1] << ' ' << (*vert)[2] << '\n';
            vert->globalNum = index++;
        }
    }

    #ifdef DEV_DEBUG
    for (auto& crys : m_crystallites)
    {
        for (auto& f_facet : crys->frontFacets())
        {
            auto facet = f_facet->facet;

            std::vector<size_t> gl_nums;
            for (auto& edge : facet->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalNum) == gl_nums.end())
                        gl_nums.push_back(vert->globalNum);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
    #else
    for (auto& sfacet : m_shellFacets)
    {
        for (auto& facet : sfacet->innerFacets())
        {
            std::vector<size_t> gl_nums;
            for (auto& edge : facet->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalNum) == gl_nums.end())
                        gl_nums.push_back(vert->globalNum);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }

    for (auto& crys : m_crystallites)
    {
        for (auto& facet : crys->innerFacets())
        {
            std::vector<size_t> gl_nums;
            for (auto& edge : facet->edges)
                for (auto& vert : edge->verts)
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalNum) == gl_nums.end())
                        gl_nums.push_back(vert->globalNum);

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
    #endif

    file.close();
}


void Polycrystal::outputLSDynaKeyword_PART(std::ofstream& file) const
{
    for (size_t id = 1; id < m_crystallites.size() + 1; id++)
    {
        file << "*PART\n"
            "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n";
        file << std::setw(10) << id;
        file << "         0";
        file << std::setw(10) << id;
        file << "         0         0         0         0         0\n";
    }
}


void Polycrystal::outputLSDynaKeyword_NODE(std::ofstream& file) const
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
    for (auto& sfacet : m_shellFacets)
    {
        for (auto& vert : sfacet->innerVerts())
        {
            vert->globalNum = index++;
            file << std::setw(8)  << vert->globalNum;
            file << std::setw(16) << (*vert)[0];
            file << std::setw(16) << (*vert)[1];
            file << std::setw(16) << (*vert)[2];
            file << "       0       0\n";
        }
    }

    for (auto& crys : m_crystallites)
    {
        for (auto& vert : crys->innerVerts())
        {
            vert->globalNum = index++;
            file << std::setw(8) << vert->globalNum;
            file << std::setw(16) << (*vert)[0];
            file << std::setw(16) << (*vert)[1];
            file << std::setw(16) << (*vert)[2];
            file << "       0       0\n";
        }
    }
}


void Polycrystal::outputLSDynaKeyword_ELEMENT_SOLID(std::ofstream& file, unsigned polycrystalId) const
{
    file << "*ELEMENT_SOLID\n"
        "$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n";
    size_t pid = 1;
    size_t eid = 10000000 * polycrystalId + 1;
    for (auto& crys : m_crystallites)
    {
        for (auto& tetr : crys->innerTetrs())
        {
            file << std::setw(8) << eid++;
            file << std::setw(8) << pid;
            file << std::setw(8) << tetr->verts[0]->globalNum;
            tva::Vec v0 = tetr->verts[1]->pos() - tetr->verts[0]->pos();
            tva::Vec v1 = tetr->verts[2]->pos() - tetr->verts[0]->pos();
            tva::Vec v2 = tetr->verts[3]->pos() - tetr->verts[0]->pos();
            if (tva::Vec::dot(v2, tva::Vec::cross(v0, v1)) > 0.0)
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


void Polycrystal::outputLSDynaKeyword(const std::string& filename, unsigned polycrystalId) const
{
    std::ofstream file(filename);
    file.setf(std::ios::right | std::ios::fixed);

    outputLSDynaKeyword_PART(file);
    outputLSDynaKeyword_NODE(file);
    outputLSDynaKeyword_ELEMENT_SOLID(file, polycrystalId);
    file << "*END";

    file.close();
}




std::string Polycrystal::generateMesh_generateLogFileName(const std::string& logFileName) const
{
    if (logFileName != "_AUTO_")
        return logFileName;

    std::stringstream ss;
    ss << "pmg_log_" << m_crystallites.size() << "_nc_";
    size_t av_nfe = 0;
    for (auto& crys : m_crystallites)
        av_nfe += crys->innerTetrs().size();
    av_nfe /= m_crystallites.size();
    ss << av_nfe << "_cfe.log";

    return ss.str();
}


void Polycrystal::generateMesh(double preferredLength, const std::string& logFileName)
{
    m_preferredLength = preferredLength;

    tva::Timer tmr;
    tmr.start();
    try
    {
        triangulateShell();
    }
    catch (std::logic_error error)
    {
        output(FileType::Obj, "debug.obj");
        throw error;
    }
    tmr.stop();

    tmr.start();
    double min_q = 0.0, av_q = 0.0;
    size_t n_elems = 0;
    size_t max = m_crystallites.size();
    #pragma omp parallel for
    for (size_t i = 0; i < max; i++)
    {
        try
        {
            m_crystallites[i]->generateMesh(m_preferredLength);
        }
        catch (std::logic_error error)
        {
            output(FileType::Obj, "debug.obj");
            throw error;
        }

        auto buf_min_av_q = m_crystallites[i]->analyzeMeshQuality();
        #pragma omp critical
        {
            min_q += buf_min_av_q.first;
            av_q += buf_min_av_q.second;
            n_elems += m_crystallites[i]->innerTetrs().size();
        }
    }
    min_q /= m_crystallites.size();
    av_q /= m_crystallites.size();
    tmr.stop();
    m_lastLogger.reset(new tva::Logger(generateMesh_generateLogFileName(logFileName)));
    *m_lastLogger << std::fixed
        << "Minimum quality"          << min_q << ""
        << "Average quality"          << av_q  << ""
        << "Crystallites number"      << m_crystallites.size() << ""
        << "Elements number"          << n_elems               << ""
        << "Preferred edge length"    << m_preferredLength     << ""
        << "Shell triangulation time" << tmr.getDuration(0,  tva::Timer::TimeScale::Microseconds) * 1e-6 << "s"
        << "Volume exhaustion time"   << tmr.getLastDuration(tva::Timer::TimeScale::Microseconds) * 1e-6 << "s";
}


const PolyMesh* Polycrystal::structurizeMesh()
{
    m_lastMesh = new PolyMesh;

    m_lastMesh->nCryses = m_crystallites.size();
    m_lastMesh->nCrysesTetrs = new size_t[m_lastMesh->nCryses];

    size_t nodes_num = m_shellVerts.size();
    for (auto& sedge : m_shellEdges)
        nodes_num += sedge->innerVerts().size();
    for (auto& sfacet : m_shellFacets)
        nodes_num += sfacet->innerVerts().size();

    for (auto& crys : m_crystallites)
        nodes_num += static_cast<size_t>(std::count_if(
            crys->innerVerts().begin(),
            crys->innerVerts().end(),
            [](auto vert) { return static_cast<bool>(vert); }));

    m_lastMesh->nNodes = nodes_num;
    m_lastMesh->nodesPositions = new double[3 * m_lastMesh->nNodes];
    Vertex** verts_ptrs = new Vertex*[m_lastMesh->nNodes];

    auto lastMesh_buf = m_lastMesh;
    auto FindNodeIndex = [lastMesh_buf, verts_ptrs](Vertex* node_ptr)-> size_t
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
    for (auto& sfacet : m_shellFacets)
    {
        for (auto& vert : sfacet->innerVerts())
        {
            verts_ptrs[vert_ind] = vert;
            m_lastMesh->nodesPositions[3 * vert_ind]     = (*vert)[0];
            m_lastMesh->nodesPositions[3 * vert_ind + 1] = (*vert)[1];
            m_lastMesh->nodesPositions[3 * vert_ind + 2] = (*vert)[2];
            vert_ind++;
        }
    }
    for (auto& crys : m_crystallites)
    {
        for (auto& vert : crys->innerVerts())
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
            m_crystallites[i]->innerTetrs().begin(),
            m_crystallites[i]->innerTetrs().end(),
            [](auto tetr) { return static_cast<bool>(tetr); }));

        tetrs_num += m_lastMesh->nCrysesTetrs[i];
    }

    m_lastMesh->nTetrs = tetrs_num;
    m_lastMesh->tetrs = new size_t[4 * m_lastMesh->nTetrs];

    size_t tetr_ind = 0;
    for (size_t i = 0; i < m_lastMesh->nCryses; i++)
    {
        for (auto& tetr : m_crystallites[i]->innerTetrs())
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


const PolyMesh* Polycrystal::getLastMesh()
{
    return m_lastMesh;
}




void Polycrystal::input(const std::string& polyStructFileName)
{
    std::ifstream input(polyStructFileName);

    size_t nodes_num;
    input >> nodes_num;

    for (size_t i = 0; i < nodes_num; i++)
    {
        double coors[3];
        input >> coors[0];
        input >> coors[1];
        input >> coors[2];

        m_shellVerts.push_back(new shell::Vertex(coors[0], coors[1], coors[2]));
    }

    size_t facets_num;
    input >> facets_num;

    for (size_t i = 0; i < facets_num; i++)
    {
        size_t facet_nodes_inds[3];
        input >> facet_nodes_inds[0];
        input >> facet_nodes_inds[1];
        input >> facet_nodes_inds[2];

        if (!findShellEdge(m_shellVerts[facet_nodes_inds[0]], m_shellVerts[facet_nodes_inds[1]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[facet_nodes_inds[0]], m_shellVerts[facet_nodes_inds[1]]));
        }

        if (!findShellEdge(m_shellVerts[facet_nodes_inds[1]], m_shellVerts[facet_nodes_inds[2]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[facet_nodes_inds[1]], m_shellVerts[facet_nodes_inds[2]]));
        }

        if (!findShellEdge(m_shellVerts[facet_nodes_inds[2]], m_shellVerts[facet_nodes_inds[0]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[facet_nodes_inds[2]], m_shellVerts[facet_nodes_inds[0]]));
        }

        m_shellFacets.push_back(new shell::Facet(
            findShellEdge(m_shellVerts[facet_nodes_inds[0]], m_shellVerts[facet_nodes_inds[1]]),
            findShellEdge(m_shellVerts[facet_nodes_inds[1]], m_shellVerts[facet_nodes_inds[2]]),
            findShellEdge(m_shellVerts[facet_nodes_inds[2]], m_shellVerts[facet_nodes_inds[0]])));
    }

    size_t cryses_num;
    input >> cryses_num;
    m_crystallites.insert(m_crystallites.end(), cryses_num, new Crystallite);

    size_t* cryses_facets_nums = new size_t[cryses_num];
    for (size_t i = 0; i < cryses_num; i++)
        input >> cryses_facets_nums[i];

    for (size_t i = 0; i < cryses_num; i++)
    {
        for (size_t j = 0; j < cryses_facets_nums[i]; j++)
        {
            size_t facet_ind;
            input >> facet_ind;

            for (int k = 0; k < 3; ++k)
            {
                if (!m_crystallites[i]->shellContains(m_shellFacets[facet_ind]->edges[k]))
                {
                    m_crystallites[i]->addToShell(m_shellFacets[facet_ind]->edges[k]);

                    if (!m_crystallites[i]->shellContains(m_shellFacets[facet_ind]->edges[k]->verts[0]))
                        m_crystallites[i]->addToShell(m_shellFacets[facet_ind]->edges[k]->verts[0]);

                    if (!m_crystallites[i]->shellContains(m_shellFacets[facet_ind]->edges[k]->verts[1]))
                        m_crystallites[i]->addToShell(m_shellFacets[facet_ind]->edges[k]->verts[1]);
                }
            }

            m_crystallites[i]->addToShell(m_shellFacets[facet_ind]);
        }
    }

    input.close();
    delete[] cryses_facets_nums;
}


void Polycrystal::input(const polygen::PolyStruct& polyStruct)
{
    for (size_t i = 0; i < polyStruct.nodes.size(); i++)
    {
        double coors[3];
        coors[0] = polyStruct.nodes[i][0];
        coors[1] = polyStruct.nodes[i][1];
        coors[2] = polyStruct.nodes[i][2];

        m_shellVerts.push_back(new shell::Vertex(coors[0], coors[1], coors[2]));
    }

    for (const auto& facet : polyStruct.facets)
    {
        if (!findShellEdge(m_shellVerts[facet[0]], m_shellVerts[facet[1]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[facet[0]], m_shellVerts[facet[1]]));
        }

        if (!findShellEdge(m_shellVerts[facet[1]], m_shellVerts[facet[2]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[facet[1]], m_shellVerts[facet[2]]));
        }

        if (!findShellEdge(m_shellVerts[facet[2]], m_shellVerts[facet[0]]))
        {
            m_shellEdges.push_back(new shell::Edge(m_shellVerts[facet[2]], m_shellVerts[facet[0]]));
        }

        m_shellFacets.push_back(new shell::Facet(
            findShellEdge(m_shellVerts[facet[0]], m_shellVerts[facet[1]]),
            findShellEdge(m_shellVerts[facet[1]], m_shellVerts[facet[2]]),
            findShellEdge(m_shellVerts[facet[2]], m_shellVerts[facet[0]])));
    }

    m_crystallites.reserve(polyStruct.cryses.size());
    for (size_t i = 0; i < polyStruct.cryses.size(); ++i)
        m_crystallites.push_back(new Crystallite(this));

    for (size_t i = 0; i < polyStruct.cryses.size(); ++i)
    {
        for (size_t j = 0; j < polyStruct.cryses[i].size(); ++j)
        {
            size_t facet_ind = polyStruct.cryses[i][j];

            for (int k = 0; k < 3; ++k)
            {
                if (!m_crystallites[i]->shellContains(m_shellFacets[facet_ind]->edges[k]))
                {
                    m_crystallites[i]->addToShell(m_shellFacets[facet_ind]->edges[k]);

                    if (!m_crystallites[i]->shellContains(m_shellFacets[facet_ind]->edges[k]->verts[0]))
                        m_crystallites[i]->addToShell(m_shellFacets[facet_ind]->edges[k]->verts[0]);

                    if (!m_crystallites[i]->shellContains(m_shellFacets[facet_ind]->edges[k]->verts[1]))
                        m_crystallites[i]->addToShell(m_shellFacets[facet_ind]->edges[k]->verts[1]);
                }
            }

            m_crystallites[i]->addToShell(m_shellFacets[facet_ind]);
        }
    }
}


std::string Polycrystal::output_generateFilename(FileType filetype, const std::string& filename) const
{
    if (filename != "_AUTO_")
        return filename;

    std::stringstream ss;
    ss << "plcr_" << m_crystallites.size() << "_nc_";
    size_t av_nfe = 0;
    for (auto& crys : m_crystallites)
        av_nfe += crys->innerTetrs().size();
    av_nfe /= m_crystallites.size();
    ss << av_nfe << "_cfe";
    switch (filetype)
    {
    case FileType::Obj:
        return ss.str() + ".obj";

    case FileType::LsDynaKeyword:
        return ss.str() + ".kw";
    }
    throw std::exception();
}


void Polycrystal::output(FileType filetype, const std::string& filename, unsigned polycrystalId) const
{
    tva::Timer tmr;
    tmr.start();
    switch (filetype)
    {
    case FileType::Obj:
        outputObj(output_generateFilename(FileType::Obj, filename));
        break;

    case FileType::LsDynaKeyword:
        outputLSDynaKeyword(output_generateFilename(FileType::LsDynaKeyword, filename), polycrystalId);
        break;
    }
    tmr.stop();

    if (m_lastLogger && m_lastLogger->isOpen())
    {
        *m_lastLogger << "Mesh file writing time" << tmr.getLastDuration(tva::Timer::TimeScale::Microseconds) * 1e-6 << "s";
        m_lastLogger->close();
    }
}




Polycrystal::Polycrystal() {}


Polycrystal::Polycrystal(const std::string& polyStructFileName)
{
    input(polyStructFileName);
}


Polycrystal::Polycrystal(const polygen::PolyStruct& polyStruct)
{
    input(polyStruct);
}


Polycrystal::~Polycrystal()
{
    for (auto& crys : m_crystallites)
        delete crys;

    for (auto& facet : m_shellFacets)
        delete facet;
    for (auto& edge : m_shellEdges)
        delete edge;
    for (auto& vert : m_shellVerts)
        delete vert;
}
