#include <cstddef>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <array>
#include <vector>
#include <map>
#include <tuple>
#include <optional>
#include "polyspt/polyhedron.h"
#include "polyspt/simplex.h"
#include "polyspt/mesh-simplex.h"
#include "mesher/algs.h"

using real_t = spt::vec<3>::real_type;
using nid_t = std::size_t;
using eid_t = std::size_t;
using pid_t = std::size_t;
using node_t = std::pair<nid_t, spt::vertex<>*>;
using element_solid_t = std::tuple<eid_t, pid_t, std::array<nid_t, 4>>;
using mesh_t = spt::mesh_v<spt::polyhedron<>, spt::elem_shape::simplex>;

enum class sections
{
    node,
    element_solid
};


std::string section_head(sections section)
{
    switch (section)
    {
    case sections::node:          return "*NODE";
    case sections::element_solid: return "*ELEMENT_SOLID";
    default: return "";
    }
}

template <typename SectionType>
SectionType parse_line(const std::string& line);

template <>
node_t parse_line(const std::string& line)
{
    std::string local_line = line;
    std::size_t pos_in_line;
    nid_t nid = std::stoull(local_line, &pos_in_line);
    spt::vec<3> pos;
    for (std::size_t i = 0; i < 3; i++)
        pos.x[i] = std::stof(local_line = local_line.substr(pos_in_line), &pos_in_line);
    return { nid, new spt::vertex<>(pos) };
}

template <>
element_solid_t parse_line(const std::string& line)
{
    std::string local_line = line;
    std::size_t pos_in_line;
    eid_t eid = std::stoull(local_line, &pos_in_line);
    pid_t pid = std::stoull(local_line = local_line.substr(pos_in_line), &pos_in_line);
    std::array<nid_t, 4> nids;
    for (std::size_t i = 0; i < 4; i++)
        nids[i] = std::stoull(local_line = local_line.substr(pos_in_line), &pos_in_line);
    return { eid, pid, nids };
}

template <typename SectionType>
SectionType parse_line(std::istream& stream)
{
    std::string line;
    std::getline(stream, line);
    return parse_line<SectionType>(line);
}

std::map<nid_t, spt::vertex<>*> parse_node_section(std::istream& stream)
{
    std::map<nid_t, spt::vertex<>*> res;
    while (!stream.eof())
    {
        if (stream.peek() == '*')
            break;
        if (stream.peek() == '$')
        {
            std::string line;
            std::getline(stream, line);
            continue;
        }
        res.insert(parse_line<node_t>(stream));
    }

    return res;
}

std::vector<element_solid_t> parse_element_solid_section(std::istream& stream)
{
    std::vector<element_solid_t> res;
    while (!stream.eof())
    {
        if (stream.peek() == '*') 
            break;
        if (stream.peek() == '$')
        {
            std::string line;
            std::getline(stream, line);
            continue;
        }
        res.push_back(parse_line<element_solid_t>(stream));
    }

    return res;
}

template <std::size_t N>
std::optional<sections> find_any_section_of(std::istream& stream, const std::array<sections, N>& sects)
{
    std::string line;
    while (std::getline(stream, line))
        if (line.front() == '*')
            for (const auto& section : sects)
                if (line.find(section_head(section)) != std::string::npos)
                    return section;

    return std::nullopt;
}

mesh_t simplices_from_kw(std::istream & stream)
{
    std::map<nid_t, spt::vertex<>*> nodes;
    std::vector<element_solid_t> elements_solid;
    while (!stream.eof())
    {
        auto section_found = find_any_section_of(stream, std::array{ sections::node, sections::element_solid });
        if (!section_found)
            break;

        switch (section_found.value())
        {
        case sections::node:
            nodes = parse_node_section(stream);
            break;
        case sections::element_solid:
            elements_solid = parse_element_solid_section(stream);
            break;
        }
    }

    std::vector<spt::vertex<>*> vertices;
    vertices.reserve(nodes.size());
    for (const auto& [key, val] : nodes)
        vertices.push_back(val);

    std::vector<spt::simplex_v<3>*> simplices;
    simplices.reserve(elements_solid.size());
    for (const auto& elem_solid : elements_solid)
        simplices.push_back(new spt::simplex_v<3>(
            nodes[std::get<2>(elem_solid)[0]],
            nodes[std::get<2>(elem_solid)[1]],
            nodes[std::get<2>(elem_solid)[2]],
            nodes[std::get<2>(elem_solid)[3]]));

    mesh_t res;
    res.vertices = std::move(vertices);
    res.elements = std::move(simplices);
    return res;
}

// Returns minimum and average quality.
std::pair<real_t, real_t> quality(const mesh_t& mesh)
{
    std::vector<real_t> qualities(mesh.elements.size());
    std::transform(mesh.elements.begin(), mesh.elements.end(), qualities.begin(), pmg::quality<3, real_t>);

    real_t min_q = *std::min_element(qualities.begin(), qualities.end());
    real_t sum_q = std::accumulate(qualities.begin(), qualities.end(), static_cast<real_t>(0.0));
    real_t av_q = sum_q / qualities.size();
    return { min_q, av_q };
}


int main()
{
    std::ifstream kw_file("heat.k");

    auto mesh = simplices_from_kw(kw_file);
    auto quality_min_av = quality(mesh);

    std::cout
        << ">>  Minimum quality....>>  " << quality_min_av.first << std::endl
        << ">>  Average quality....>>  " << quality_min_av.second << std::endl
        << ">>  Elements number....>>  " << mesh.elements.size() << std::endl;

    return 0;
}
