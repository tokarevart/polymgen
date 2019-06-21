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
#include "polyspt/simplex.h"

using real_t = spt::vec<3>::real_type;
using nid_t = std::size_t;
using eid_t = std::size_t;
using pid_t = std::size_t;
using node_t = std::pair<nid_t, spt::vertex<>*>;
using element_solid_t = std::tuple<eid_t, pid_t, std::array<nid_t, 4>>;

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
    std::size_t pos_in_line;
    nid_t nid = std::stoull(line, &pos_in_line);
    spt::vec<3> pos;
    for (std::size_t i = 0; i < 3; i++)
        pos.x[i] = std::stof(line.substr(pos_in_line), &pos_in_line);
    return { nid, new spt::vertex<>(pos) };
}

template <>
element_solid_t parse_line(const std::string& line)
{
    std::size_t pos_in_line;
    eid_t eid = std::stoull(line, &pos_in_line);
    pid_t pid = std::stoull(line, &pos_in_line);
    std::array<nid_t, 4> nids;
    for (std::size_t i = 0; i < 4; i++)
        nids[i] = std::stoull(line.substr(pos_in_line), &pos_in_line);
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

std::vector<spt::simplex_v<3>*> simplices_from_kw(std::istream& stream)
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

    std::vector<spt::simplex_v<3>*> simplices;
    simplices.reserve(elements_solid.size());
    for (const auto& elem_solid : elements_solid)
        simplices.push_back(new spt::simplex_v<3>(
            nodes[std::get<2>(elem_solid)[0]],
            nodes[std::get<2>(elem_solid)[1]],
            nodes[std::get<2>(elem_solid)[2]],
            nodes[std::get<2>(elem_solid)[3]]));

    return simplices;
}


int main()
{
    std::ifstream kw_file("mesh.kw");
    auto simplices = simplices_from_kw(kw_file);

    return 0;
}
