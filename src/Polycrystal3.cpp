#include "Polycrystal3.h"
#include <algorithm>
#include <sstream>
#include <iomanip>
#include <iostream>
#include "helpers/spatialalgs/SpatialAlgs.h"
#include "helpers/Timer.h"

//#define DEV_DEBUG




ShellEdge3* Polycrystal3::findShellEdge(const ShellVertex3* v0, const ShellVertex3* v1) const
{
    for (auto& s_edge : _shellEdges)
    {
        if ((s_edge->verts[0] == v0  &&
             s_edge->verts[1] == v1) ||
            (s_edge->verts[1] == v0  &&
             s_edge->verts[0] == v1))
            return s_edge;
    }
    
    return nullptr;
}

Edge3* Polycrystal3::findStartFrontEdge(const Vertex3* v0, const Vertex3* v1) const
{
    for (auto& edge : _startFrontEdges)
    {
        if ((edge->verts[0] == v0  &&
             edge->verts[1] == v1) ||
            (edge->verts[1] == v0  &&
             edge->verts[0] == v1))
            return edge;
    }

    return nullptr;
}




void Polycrystal3::triangulateShell()
{
    for (size_t i = 0ull; i < _startFrontEdges.size();)
    {
        Edge3* cur_edge = *std::next(_startFrontEdges.begin(), i);
        if (cur_edge->sqrMagnitude() > 2.25 * m_preferredLength * m_preferredLength)
        {
            cur_edge->make2Instead(_startFrontFacets, _startFrontEdges, _startFrontVertexes);
        }
        else i++;
    }
}


void Polycrystal3::setLinksWithShell()
{
    double sqr_sufficient_dist = m_preferredLength * m_preferredLength * 1e-6;
    size_t shell_verts_num = _shellVertexes.size();
    size_t shell_edges_num = _shellEdges.size();
    size_t shell_facets_num = _shellFacets.size();
    #pragma omp parallel firstprivate(sqr_sufficient_dist, shell_verts_num, shell_edges_num, shell_facets_num)
    {
    #pragma omp for
        for (size_t i = 0ull; i < shell_verts_num; i++)
        {
            for (auto& vert : _startFrontVertexes)
            {
                if ((*_shellVertexes[i] - *vert).sqrMagnitude() < sqr_sufficient_dist)
                {
                    vert->belongsToShellVertex = _shellVertexes[i];
                    break;
                }
            }
        }
        #pragma omp for
        for (size_t i = 0ull; i < shell_edges_num; i++)
        {
            tva::Point3 proj_buf;
            for (auto& vert : _startFrontVertexes)
            {
                if (vert->belongsToShellVertex)
                    continue;

                if (tva::spatialalgs::project(
                    proj_buf,
                    vert->getPos(),
                    _shellEdges[i]->verts[0]->getPos(), _shellEdges[i]->verts[1]->getPos()) &&
                    (proj_buf - vert->getPos()).sqrMagnitude() < sqr_sufficient_dist)
                {
                    vert->belongsToShellEdge = _shellEdges[i];
                }
            }
        }
        #pragma omp for
        for (size_t i = 0ull; i < shell_facets_num; i++)
        {
            for (auto& vert : _startFrontVertexes)
            {
                if (vert->belongsToShellVertex ||
                    vert->belongsToShellEdge)
                    continue;

                tva::Vec3 v_pos = vert->getPos();
                tva::Vec3 v0_pos = _shellFacets[i]->edges[0]->verts[0]->getPos();
                tva::Vec3 v1_pos = _shellFacets[i]->edges[0]->verts[1]->getPos();
                tva::Vec3 third_pos = _shellFacets[i]->findVertNotIncludedInEdge(_shellFacets[i]->edges[0])->getPos();

                if (!tva::spatialalgs::isPointOnTriangle(v_pos, v0_pos, v1_pos, third_pos))
                    continue;

                vert->belongsToShellFacet = _shellFacets[i];
            }
        }
    }
}


void Polycrystal3::startFrontDelaunayPostprocessing()
{
    Vertex3* around_nodes[2];
    Facet3* around_facets[2];
    size_t edges_num = _startFrontEdges.size();
    for (size_t i = 0ull; i < edges_num; i++)
    {
        Edge3* cur_edge = *std::next(_startFrontEdges.begin(), i);
        if ((cur_edge->verts[0]->belongsToShellEdge == cur_edge->verts[1]->belongsToShellEdge &&
             cur_edge->verts[1]->belongsToShellEdge) ||
            (cur_edge->verts[0]->belongsToShellEdge &&
             cur_edge->verts[1]->belongsToShellVertex &&
             cur_edge->verts[0]->belongsToShellEdge->contains(cur_edge->verts[1]->belongsToShellVertex)) ||
            (cur_edge->verts[1]->belongsToShellEdge &&
             cur_edge->verts[0]->belongsToShellVertex &&
             cur_edge->verts[1]->belongsToShellEdge->contains(cur_edge->verts[0]->belongsToShellVertex)))
            continue;
        

        cur_edge->find2AdjFacets(_startFrontFacets, around_facets[0], around_facets[1]);
        around_nodes[0] = around_facets[0]->findVertNotIncludedInEdge(cur_edge);
        around_nodes[1] = around_facets[1]->findVertNotIncludedInEdge(cur_edge);
        if (around_nodes[0]->belongsToShellEdge && around_nodes[1]->belongsToShellVertex ||
            around_nodes[1]->belongsToShellEdge && around_nodes[0]->belongsToShellVertex)
            continue;

        if (cur_edge->flipIfNeeded(_startFrontEdges, _startFrontFacets))
            i--;
    }
}




void Polycrystal3::outputDataObj(const std::string& filename) const
{
    std::ofstream file(filename);

    size_t index = 1ull;
    for (auto& vert : _startFrontVertexes)
    {

        file << "v " << (*vert)[0] << ' ' << (*vert)[1] << ' ' << (*vert)[2] << '\n';
        vert->globalNum = index++;
    }
    for (auto& crys : _crystallites)
    {
        if (!crys)
            continue;

        for (auto& vert : crys->getInnerVertexes())
        {
            file << "v " << (*vert)[0] << ' ' << (*vert)[1] << ' ' << (*vert)[2] << '\n';
            vert->globalNum = index++;
        }
    }

#ifdef DEV_DEBUG
    for (auto& crys : _crystallites)
    {
        for (auto& f_facet : crys->m_frontFacets)
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
    for (auto& facet : _startFrontFacets)
    {
        std::vector<size_t> gl_nums;
        for (auto& edge : facet->edges)
        {
            for (auto& vert : edge->verts)
            {
                if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalNum) == gl_nums.end())
                    gl_nums.push_back(vert->globalNum);
            }
        }

        file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
    }

    for (auto& crys : _crystallites)
    {
        for (auto& facet : crys->getInnerFacets())
        {
            std::vector<size_t> gl_nums;
            for (auto& edge : facet->edges)
            {
                for (auto& vert : edge->verts)
                {
                    if (std::find(gl_nums.begin(), gl_nums.end(), vert->globalNum) == gl_nums.end())
                        gl_nums.push_back(vert->globalNum);
                }
            }

            file << "f " << gl_nums[0] << ' ' << gl_nums[1] << ' ' << gl_nums[2] << '\n';
        }
    }
#endif

    file.close();
}


void Polycrystal3::outputDataLSDynaKeyword_PART(std::ofstream& file) const
{
    for (size_t id = 1ull, max_id = _crystallites.size() + 1ull; id < max_id; id++)
    {
        file << "*PART\n"
            "$#     pid     secid       mid     eosid      hgid      grav    adpopt      tmid\n";
        file << std::setw(10) << id;
        file << "         0";
        file << std::setw(10) << id;
        file << "         0         0         0         0         0\n";
    }
}


void Polycrystal3::outputDataLSDynaKeyword_NODE(std::ofstream& file) const
{
    file << "*NODE\n"
        "$#   nid               x               y               z      tc      rc\n";
    size_t index = 1ull;
    for (auto& vert : _startFrontVertexes)
    {
        vert->globalNum = index++;
        file << std::setw(8)  << vert->globalNum;
        file << std::setw(16) << (*vert)[0];
        file << std::setw(16) << (*vert)[1];
        file << std::setw(16) << (*vert)[2];
        file << "       0       0\n";
    }
    for (auto& crys : _crystallites)
    {
        for (auto& vert : crys->getInnerVertexes())
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


void Polycrystal3::outputDataLSDynaKeyword_ELEMENT_SOLID(std::ofstream& file, int polycrystalId) const
{
    file << "*ELEMENT_SOLID\n"
        "$#   eid     pid      n1      n2      n3      n4      n5      n6      n7      n8\n";
    size_t pid = 1ull;
    size_t eid = 10000000ull * polycrystalId + 1ull;
    for (auto& crys : _crystallites)
    {
        for (auto& simp : crys->getInnerSimplexes3())
        {
            file << std::setw(8) << eid++;
            file << std::setw(8) << pid;
            file << std::setw(8) << simp->verts[0]->globalNum;
            tva::Vec3 v0 = simp->verts[1]->getPos() - simp->verts[0]->getPos();
            tva::Vec3 v1 = simp->verts[2]->getPos() - simp->verts[0]->getPos();
            tva::Vec3 v2 = simp->verts[3]->getPos() - simp->verts[0]->getPos();
            if (tva::Vec3::dotProduct(v2, tva::Vec3::crossProduct(v0, v1)) > 0.0)
            {
                file << std::setw(8) << simp->verts[1]->globalNum;
                file << std::setw(8) << simp->verts[2]->globalNum;
            }
            else
            {
                file << std::setw(8) << simp->verts[2]->globalNum;
                file << std::setw(8) << simp->verts[1]->globalNum;
            }
            file << std::setw(8) << simp->verts[3]->globalNum;
            file << "       0       0       0       0\n";
        }
        pid++;
    }
}


void Polycrystal3::outputDataLSDynaKeyword(const std::string& filename, int polycrystalId) const
{
    std::ofstream file(filename);
    file.setf(std::ios::right | std::ios::fixed);

    outputDataLSDynaKeyword_PART(file);
    outputDataLSDynaKeyword_NODE(file);
    outputDataLSDynaKeyword_ELEMENT_SOLID(file);
    file << "*END";

    file.close();
}




std::string Polycrystal3::generateMesh_generateLogFileName(const std::string& logFileName) const
{
    if (logFileName != "_AUTO_")
        return logFileName;

    std::stringstream ss;
    ss << "log_" << _crystallites.size() << "_ncr_";
    size_t av_nfe = 0ull;
    for (auto& crys : _crystallites)
        av_nfe += crys->getInnerSimplexes3().size();
    av_nfe /= _crystallites.size();
    ss << av_nfe << "_fe_per_cr.log";

    return ss.str();
}


void Polycrystal3::generateMesh(const double preferredLength, const std::string& logFileName)
{
    m_preferredLength = preferredLength;

    tva::Timer tmr;
    tmr.start();
    triangulateShell();
    setLinksWithShell();
    startFrontDelaunayPostprocessing();
    tmr.stop();

    tmr.start();
    double min_q = 0.0, av_q = 0.0;
    size_t n_elems = 0ull;
    #pragma omp parallel for
    for (size_t i = 0ull, max = _crystallites.size(); i < max; i++)
    {
        _crystallites[i]->setStartFront(_startFrontEdges, _startFrontFacets);
        //try 
        //{ 
            _crystallites[i]->generateMesh(m_preferredLength); 
        //}
        //catch (std::logic_error error)
        //{
        //    outputData(OBJ, "debug.obj");
        //    throw error;
        //}

        auto buf_min_av_q = _crystallites[i]->analyzeMeshQuality();
        #pragma omp critical
        {
            min_q += buf_min_av_q.first;
            av_q += buf_min_av_q.second;
            n_elems += _crystallites[i]->getInnerSimplexes3().size();
        }
    }
    min_q /= _crystallites.size();
    av_q /= _crystallites.size();
    tmr.stop();
    _lastLogger.reset(new tva::Logger(generateMesh_generateLogFileName(logFileName)));
    *_lastLogger << std::fixed
        << "Minimum quality"          << min_q << ""
        << "Average quality"          << av_q  << ""
        << "Crystallites number"      << _crystallites.size() << ""
        << "Elements number"          << n_elems              << ""
        << "Shell triangulation time" << tmr.getDuration(0,  tva::Timer::MICROSECONDS) * 1e-6 << "s"
        << "Volume exhaustion time"   << tmr.getLastDuration(tva::Timer::MICROSECONDS) * 1e-6 << "s";
}


void Polycrystal3::generateMesh(const std::string& polyStructFileName, const double preferredLength, const std::string& logFileName)
{
    inputData(polyStructFileName);
    generateMesh(preferredLength, logFileName);
}


void Polycrystal3::generateMesh(const PolyStruct* polyStruct, const double preferredLength, const std::string& logFileName)
{
    inputData(polyStruct);
    generateMesh(preferredLength, logFileName);
}


const PolyMesh* Polycrystal3::structurizeMesh()
{
    _lastMesh = new PolyMesh;

    _lastMesh->nCryses = _crystallites.size();
    _lastMesh->nCrysesTetrs = new size_t[_lastMesh->nCryses];

    size_t nodes_num = _startFrontVertexes.size();
    for (auto& crys : _crystallites)
        nodes_num += std::count_if(
            crys->getInnerVertexes().begin(),
            crys->getInnerVertexes().end(),
            [](auto vert)-> bool { return (bool)vert; });

    _lastMesh->nNodes = nodes_num;
    _lastMesh->nodesPositions = new double[3ull * _lastMesh->nNodes];
    Vertex3** nodes_ptrs = new Vertex3*[_lastMesh->nNodes];

    auto lastMesh_buf = _lastMesh;
    auto FindNodeIndex = [lastMesh_buf, nodes_ptrs](Vertex3* node_ptr)-> size_t
    {
        for (size_t i = 0ull; i < lastMesh_buf->nNodes; i++)
        {
            if (nodes_ptrs[i] == node_ptr)
                return i;
        }
        return lastMesh_buf->nNodes;
    };

    size_t node_ind = 0ull;
    for (auto& s_f_vert : _startFrontVertexes)
    {
        nodes_ptrs[node_ind] = s_f_vert;
        _lastMesh->nodesPositions[3ull * node_ind]        = (*s_f_vert)[0];
        _lastMesh->nodesPositions[3ull * node_ind + 1ull] = (*s_f_vert)[1];
        _lastMesh->nodesPositions[3ull * node_ind + 2ull] = (*s_f_vert)[2];
        node_ind++;
    }
    for (auto& crys : _crystallites)
    {
        for (auto& vert : crys->getInnerVertexes())
        {
            nodes_ptrs[node_ind] = vert;
            _lastMesh->nodesPositions[3ull * node_ind]        = (*vert)[0];
            _lastMesh->nodesPositions[3ull * node_ind + 1ull] = (*vert)[1];
            _lastMesh->nodesPositions[3ull * node_ind + 2ull] = (*vert)[2];
            node_ind++;
        }
    }

    size_t tetrs_num = 0ull;
    for (size_t i = 0ull; i < _lastMesh->nCryses; i++)
    {
        _lastMesh->nCrysesTetrs[i] = std::count_if(
            _crystallites[i]->getInnerSimplexes3().begin(),
            _crystallites[i]->getInnerSimplexes3().end(),
            [](auto simp)-> bool { return (bool)simp; });

        tetrs_num += _lastMesh->nCrysesTetrs[i];
    }

    _lastMesh->nTetrs = tetrs_num;
    _lastMesh->tetrs = new size_t[4ull * _lastMesh->nTetrs];

    size_t tetr_ind = 0ull;
    for (size_t i = 0ull; i < _lastMesh->nCryses; i++)
    {
        for (auto& tetr : _crystallites[i]->getInnerSimplexes3())
        {
            _lastMesh->tetrs[4ull * tetr_ind]        = FindNodeIndex(tetr->verts[0]);
            _lastMesh->tetrs[4ull * tetr_ind + 1ull] = FindNodeIndex(tetr->verts[1]);
            _lastMesh->tetrs[4ull * tetr_ind + 2ull] = FindNodeIndex(tetr->verts[1]);
            _lastMesh->tetrs[4ull * tetr_ind + 3ull] = FindNodeIndex(tetr->verts[1]);
            tetr_ind++;
        }
    }

    return _lastMesh;
}


const PolyMesh* Polycrystal3::getLastMesh()
{
    return _lastMesh;
}




void Polycrystal3::inputData(const std::string& polyStructFileName)
{
    std::ifstream input(polyStructFileName);

    size_t nodes_num;
    input >> nodes_num;

    for (size_t i = 0ull; i < nodes_num; i++)
    {
        double coors[3];
        input >> coors[0];
        input >> coors[1];
        input >> coors[2];

        _shellVertexes.push_back(new ShellVertex3(coors[0], coors[1], coors[2]));
        _startFrontVertexes.push_back(new Vertex3(coors[0], coors[1], coors[2]));
    }

    size_t facets_num;
    input >> facets_num;

    for (size_t i = 0ull; i < facets_num; i++)
    {
        size_t facet_nodes_inds[3];
        input >> facet_nodes_inds[0];
        input >> facet_nodes_inds[1];
        input >> facet_nodes_inds[2];

        if (!findShellEdge(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]))
        {
            _shellEdges.push_back(new ShellEdge3(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]));
            _startFrontEdges.push_back(new Edge3(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[0]), 
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[1])));
        }

        if (!findShellEdge(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]))
        {
            _shellEdges.push_back(new ShellEdge3(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]));
            _startFrontEdges.push_back(new Edge3(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[1]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[2])));
        }

        if (!findShellEdge(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]]))
        {
            _shellEdges.push_back(new ShellEdge3(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]]));
            _startFrontEdges.push_back(new Edge3(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[2]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[0])));
        }

        _shellFacets.push_back(new ShellFacet3(
            findShellEdge(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]),
            findShellEdge(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]),
            findShellEdge(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]])));
        _startFrontFacets.push_back(new Facet3(
            findStartFrontEdge(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[0]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[1])),
            findStartFrontEdge(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[1]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[2])),
            findStartFrontEdge(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[2]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[0]))));
    }

    size_t cryses_num;
    input >> cryses_num;
    _crystallites.insert(_crystallites.end(), cryses_num, new Crystallite3);

    size_t* cryses_facets_nums = new size_t[cryses_num];
    for (size_t i = 0ull; i < cryses_num; i++)
        input >> cryses_facets_nums[i];

    for (size_t i = 0ull; i < cryses_num; i++)
    {
        for (size_t j = 0ull; j < cryses_facets_nums[i]; j++)
        {
            size_t facet_ind;
            input >> facet_ind;

            if (!_crystallites[i]->shellEdgesContains(_shellFacets[facet_ind]->edges[0]))
                _crystallites[i]->addShellEdge(_shellFacets[facet_ind]->edges[0]);

            if (!_crystallites[i]->shellEdgesContains(_shellFacets[facet_ind]->edges[1]))
                _crystallites[i]->addShellEdge(_shellFacets[facet_ind]->edges[1]);

            if (!_crystallites[i]->shellEdgesContains(_shellFacets[facet_ind]->edges[2]))
                _crystallites[i]->addShellEdge(_shellFacets[facet_ind]->edges[2]);

            _crystallites[i]->addShellFacet(_shellFacets[facet_ind]);
        }
    }

    input.close();
    delete[] cryses_facets_nums;
}


void Polycrystal3::inputData(const PolyStruct* polyStruct)
{
    for (size_t i = 0ull; i < polyStruct->nNodes; i++)
    {
        double coors[3];
        coors[0] = polyStruct->nodesPositions[3ull * i];
        coors[1] = polyStruct->nodesPositions[3ull * i + 1ull];
        coors[2] = polyStruct->nodesPositions[3ull * i + 2ull];

        _shellVertexes.push_back(new ShellVertex3(coors[0], coors[1], coors[2]));
        _startFrontVertexes.push_back(new Vertex3(coors[0], coors[1], coors[2]));
    }

    for (size_t i = 0ull; i < polyStruct->nFacets; i++)
    {
        size_t facet_nodes_inds[3];
        facet_nodes_inds[0] = polyStruct->facets[3ull * i];
        facet_nodes_inds[1] = polyStruct->facets[3ull * i + 1ull];
        facet_nodes_inds[2] = polyStruct->facets[3ull * i + 2ull];

        if (!findShellEdge(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]))
        {
            _shellEdges.push_back(new ShellEdge3(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]));
            _startFrontEdges.push_back(new Edge3(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[0]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[1])));
        }

        if (!findShellEdge(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]))
        {
            _shellEdges.push_back(new ShellEdge3(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]));
            _startFrontEdges.push_back(new Edge3(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[1]), 
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[2])));
        }

        if (!findShellEdge(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]]))
        {
            _shellEdges.push_back(new ShellEdge3(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]]));
            _startFrontEdges.push_back(new Edge3(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[2]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[0])));
        }

        _shellFacets.push_back(new ShellFacet3(
            findShellEdge(_shellVertexes[facet_nodes_inds[0]], _shellVertexes[facet_nodes_inds[1]]),
            findShellEdge(_shellVertexes[facet_nodes_inds[1]], _shellVertexes[facet_nodes_inds[2]]),
            findShellEdge(_shellVertexes[facet_nodes_inds[2]], _shellVertexes[facet_nodes_inds[0]])));
        _startFrontFacets.push_back(new Facet3(
            findStartFrontEdge(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[0]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[1])),
            findStartFrontEdge(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[1]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[2])),
            findStartFrontEdge(
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[2]),
                *std::next(_startFrontVertexes.begin(), facet_nodes_inds[0]))));
    }

    for (size_t i = 0ull; i < polyStruct->nCryses; i++)
        _crystallites.push_back(new Crystallite3);
    
    size_t cryses_inner_ind = 0ull;
    for (size_t i = 0ull; i < polyStruct->nCryses; i++)
    {
        for (size_t j = 0ull; j < polyStruct->nCrysesFacets[i]; j++)
        {
            size_t facet_ind = polyStruct->cryses[cryses_inner_ind++];

            if (!_crystallites[i]->shellEdgesContains(_shellFacets[facet_ind]->edges[0]))
                _crystallites[i]->addShellEdge(_shellFacets[facet_ind]->edges[0]);

            if (!_crystallites[i]->shellEdgesContains(_shellFacets[facet_ind]->edges[1]))
                _crystallites[i]->addShellEdge(_shellFacets[facet_ind]->edges[1]);

            if (!_crystallites[i]->shellEdgesContains(_shellFacets[facet_ind]->edges[2]))
                _crystallites[i]->addShellEdge(_shellFacets[facet_ind]->edges[2]);

            _crystallites[i]->addShellFacet(_shellFacets[facet_ind]);
        }
    }
}


std::string Polycrystal3::outputData_generateFilename(FileType filetype, const std::string& filename) const
{
    if (filename != "_AUTO_")
        return filename;

    std::stringstream ss;
    ss << "polycr_" << _crystallites.size() << "_ncr_";
    size_t av_nfe = 0ull;
    for (auto& crys : _crystallites)
        av_nfe += crys->getInnerSimplexes3().size();
    av_nfe /= _crystallites.size();
    ss << av_nfe << "_fe_per_cr";
    switch (filetype)
    {
    case OBJ:
        return ss.str() + ".obj";
        break;

    case LS_DYNA_KEYWORD:
        return ss.str() + ".kw";
        break;

    default:
        throw std::range_error("Wrong FileType.\nError in function: Polycrystal3::outputData_generateFilename");
    }
}


void Polycrystal3::outputData(FileType filetype, const std::string& filename, int polycrystalId) const
{
    tva::Timer tmr;
    tmr.start();
    switch (filetype)
    {
    case OBJ:
        outputDataObj(outputData_generateFilename(OBJ, filename));
        break;

    case LS_DYNA_KEYWORD:
        outputDataLSDynaKeyword(outputData_generateFilename(LS_DYNA_KEYWORD, filename), polycrystalId);
        break;

    default:
        throw std::range_error("Wrong FileType.\nError in function: Polycrystal3::outputData");
    }
    tmr.stop();

    #ifndef DEV_DEBUG
    if (_lastLogger->isOpen())
    {
        *_lastLogger << "Mesh file writing time" << tmr.getLastDuration(tva::Timer::MICROSECONDS) * 1e-6 << "s";
        _lastLogger->close();
    }
    #endif
}




Polycrystal3::Polycrystal3() {}


Polycrystal3::Polycrystal3(const std::string& polyStructFileName)
{
    inputData(polyStructFileName);
}


Polycrystal3::Polycrystal3(const PolyStruct* polyStruct)
{
    inputData(polyStruct);
}


Polycrystal3::~Polycrystal3()
{
    for (auto& crys : _crystallites)
        delete crys;

    for (auto& facet : _shellFacets)
        delete facet;
    for (auto& edge : _shellEdges)
        delete edge;
    for (auto& vert : _shellVertexes)
        delete vert;

    for (auto& f_facet : _startFrontFacets)
        delete f_facet;
    for (auto& f_edge : _startFrontEdges)
        delete f_edge;
    for (auto& f_vert : _startFrontVertexes)
        delete f_vert;
}