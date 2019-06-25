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
#include <memory>
#include "polyspt/polyhedron.h"
#include "polyspt/simplex.h"
#include "polyspt/mesh.h"
#include "mesher/algs.h"

using real_t = double;
using coordinate_t = real_t;
using vec_t = spt::vec<3, coordinate_t>;
using nid_t = std::size_t;
using eid_t = std::size_t;
using pid_t = std::size_t;
using vertex_t = spt::vertex<3, coordinate_t>;
using node_t = std::pair<nid_t, std::unique_ptr<vertex_t>>;
using tetrahedron_t = spt::tetrahedron_v<3, coordinate_t>;
using element_solid_t = std::tuple<eid_t, pid_t, std::array<nid_t, 4>>;
using mesh_t = spt::unique_mesh<tetrahedron_t>;

enum class section {
    node,
    element_solid
};


std::string section_head(section sect) {
    switch (sect) {
    case section::node:          return "*NODE";
    case section::element_solid: return "*ELEMENT_SOLID";
    default: return "";
    }
}

template <typename SectionType>
SectionType parse_line(const std::string& line);

template <>
node_t parse_line(const std::string& line) {
    std::string local_line = line;
    std::size_t pos_in_line;
    nid_t nid = std::stoull(local_line, &pos_in_line);
    vec_t pos;
    for (std::size_t i = 0; i < 3; i++)
        pos.x[i] = std::stof(local_line = local_line.substr(pos_in_line), &pos_in_line);
    return { nid, std::make_unique<vertex_t>(pos) };
}

template <>
element_solid_t parse_line(const std::string& line) {
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
SectionType parse_line(std::istream& stream) {
    std::string line;
    std::getline(stream, line);
    return parse_line<SectionType>(line);
}

std::map<nid_t, std::unique_ptr<vertex_t>> parse_node_section(std::istream& stream) {
    std::map<nid_t, std::unique_ptr<vertex_t>> res;
    while (!stream.eof()) {
        if (stream.peek() == '*')
            break;
        if (stream.peek() == '$') {
            std::string line;
            std::getline(stream, line);
            continue;
        }
        res.insert(parse_line<node_t>(stream));
    }

    return res;
}

std::vector<element_solid_t> parse_element_solid_section(std::istream& stream) {
    std::vector<element_solid_t> res;
    while (!stream.eof()) {
        if (stream.peek() == '*')
            break;
        if (stream.peek() == '$') {
            std::string line;
            std::getline(stream, line);
            continue;
        }
        res.push_back(parse_line<element_solid_t>(stream));
    }

    return res;
}

template <std::size_t N>
std::optional<section> find_any_section_of(std::istream& stream, const std::array<section, N>& sects) {
    std::string line;
    while (std::getline(stream, line))
        if (line.front() == '*')
            for (const auto& section : sects)
                if (line.find(section_head(section)) != std::string::npos)
                    return section;

    return std::nullopt;
}

mesh_t simplices_from_kw(std::istream& stream) {
    std::map<nid_t, std::unique_ptr<vertex_t>> nodes;
    std::vector<element_solid_t> elements_solid;
    while (!stream.eof()) {
        auto section_found = find_any_section_of(stream, std::array{ section::node, section::element_solid });
        if (!section_found)
            break;

        switch (section_found.value()) {
        case section::node:
            nodes = parse_node_section(stream);
            break;
        case section::element_solid:
            elements_solid = parse_element_solid_section(stream);
            break;
        }
    }

    std::vector<std::unique_ptr<tetrahedron_t>> simplices;
    simplices.reserve(elements_solid.size());
    for (const auto& elem_solid : elements_solid)
        simplices.push_back(std::make_unique<tetrahedron_t>(
            nodes[std::get<2>(elem_solid)[0]].get(),
            nodes[std::get<2>(elem_solid)[1]].get(),
            nodes[std::get<2>(elem_solid)[2]].get(),
            nodes[std::get<2>(elem_solid)[3]].get()));

    std::vector<std::unique_ptr<vertex_t>> vertices;
    vertices.reserve(nodes.size());
    for (auto&& [key, val] : nodes)
        vertices.push_back(std::move(val));

    mesh_t res;
    res.vertices = std::move(vertices);
    res.elements = std::move(simplices);
    return res;
}

// Returns minimum and average quality.
std::pair<real_t, real_t> quality(const mesh_t& mesh) {
    std::vector<real_t> qualities(mesh.elements.size());
    std::transform(mesh.elements.begin(), mesh.elements.end(), qualities.begin(), 
                   [](const auto& elem) { return pmg::quality(elem.get()); });

    real_t min_q = *std::min_element(qualities.begin(), qualities.end());
    real_t sum_q = std::accumulate(qualities.begin(), qualities.end(), static_cast<real_t>(0.0));
    real_t av_q = sum_q / qualities.size();
    return { min_q, av_q };
}

template <template <typename Arg> typename T, typename Arg>
struct A;

template <template <typename Arg> typename T>
struct A<T, int> {
    T<int> a;
};

template <typename T>
using raw_ptr = T*;

template <typename Arg>
using Araw = A<raw_ptr, Arg>;

int main() {
    std::ifstream kw_file("heat.k");

    auto mesh = simplices_from_kw(kw_file);
    auto quality_min_av = quality(mesh);

    std::cout
        << ">>  Minimum quality....>>  " << quality_min_av.first << std::endl
        << ">>  Average quality....>>  " << quality_min_av.second << std::endl
        << ">>  Elements number....>>  " << mesh.elements.size() << std::endl;
        
    std::cout << spt::dot(spt::vec<3>(0, 1, 1), spt::vec<3>(0, 1, 2)) << std::endl;
    std::cout << spt::dot(spt::mat<3>::identity(), spt::vec<3>(0, 1, 2)).magnitude() << std::endl;
    std::cout << spt::dot(spt::mat<3>::identity().inversed(), spt::mat<3>::identity().transposed() * 2)[1].magnitude() << std::endl;
    A<raw_ptr, int> kek;
    kek.a = new int(42);
    std::cout << *kek.a << std::endl;

    return 0;
}
