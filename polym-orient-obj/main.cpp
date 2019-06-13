#include <cstddef>
#include <iostream>
#include <algorithm>
#include <array>
#include <string>
#include <chrono>
#include <optional>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include "vec.h"


using indices = std::vector<tinyobj::index_t>;


bool faceContainsVert( indices::const_iterator face, indices::const_iterator vert )
{
    if (   (face + 0)->vertex_idx == vert->vertex_idx
        || (face + 1)->vertex_idx == vert->vertex_idx
        || (face + 2)->vertex_idx == vert->vertex_idx)
        return true;

    return false;
}


bool findOneVertNot(std::size_t& out, indices::const_iterator in_face, indices::const_iterator not_face)
{
    std::size_t n_not_conts = 0;
    std::array<bool, 3> not_conts;
    if (not_conts[0] = !faceContainsVert(not_face, in_face + 0)) n_not_conts++;
    if (not_conts[1] = !faceContainsVert(not_face, in_face + 1)) n_not_conts++;
    if (not_conts[2] = !faceContainsVert(not_face, in_face + 2)) n_not_conts++;
    if (n_not_conts != 1) return false;

    if (not_conts[0])
    {
        out = (in_face + 0)->vertex_idx;
        return true;
    }
    if (not_conts[1])
    {
        out = (in_face + 1)->vertex_idx;
        return true;
    }
    if (not_conts[2])
    {
        out = (in_face + 2)->vertex_idx;
        return true;
    }

    throw std::logic_error("Function findOneVertNot didn't find vertex.");
}


std::optional<std::size_t> opp_vert(indices::const_iterator face, const tinyobj::shape_t& shape)
{
    std::vector<std::size_t> opp_idces;

    for (std::size_t i = 0; i < shape.mesh.indices.size(); i += 3)
    {
        std::size_t opp_idx;
        if (findOneVertNot(opp_idx, shape.mesh.indices.begin() + i, face))
        {
            if (std::count(opp_idces.begin(), opp_idces.end(), opp_idx) == 2)
                return std::optional(opp_idx);

            opp_idces.push_back(opp_idx);
        }
    }

    return std::nullopt;
}


spt::vec3 vert_pos(std::size_t idx, const tinyobj::attrib_t& attrib)
{
    return spt::vec3(attrib.vertices[3 * idx + 0],
                     attrib.vertices[3 * idx + 1],
                     attrib.vertices[3 * idx + 2]);
}


spt::vec3 raw_normal(indices::const_iterator face, const tinyobj::attrib_t& attrib)
{
    spt::vec3 v0(attrib.vertices[3 * (face + 0)->vertex_idx + 0],
                 attrib.vertices[3 * (face + 0)->vertex_idx + 1],
                 attrib.vertices[3 * (face + 0)->vertex_idx + 2]);

    spt::vec3 v1(attrib.vertices[3 * (face + 1)->vertex_idx + 0],
                 attrib.vertices[3 * (face + 1)->vertex_idx + 1],
                 attrib.vertices[3 * (face + 1)->vertex_idx + 2]);

    spt::vec3 v2(attrib.vertices[3 * (face + 2)->vertex_idx + 0],
                 attrib.vertices[3 * (face + 2)->vertex_idx + 1],
                 attrib.vertices[3 * (face + 2)->vertex_idx + 2]);

    return spt::vec3::cross(v1 - v0, v2 - v0);
}


void orientate_mesh(const tinyobj::attrib_t& attrib, tinyobj::shape_t& shape, std::vector<size_t>& dont_write_faces)
{
    std::vector<std::pair<indices::iterator, indices::iterator>> indices_to_swap;

    size_t n_indices = shape.mesh.indices.size();
    #pragma omp parallel for
    for (std::size_t idx_offset = 0; idx_offset < n_indices; idx_offset += 3)
    {
        auto opp_vert_idx = opp_vert(shape.mesh.indices.begin() + idx_offset, shape);
        if (!opp_vert_idx)
            #pragma omp critical
            dont_write_faces.push_back(idx_offset);
        spt::vec3 opp_vert_pos = vert_pos(*opp_vert_idx, attrib);

        spt::vec3 v0_to_opp = opp_vert_pos - spt::vec3(attrib.vertices[3 * shape.mesh.indices[idx_offset].vertex_idx + 0],
                                                       attrib.vertices[3 * shape.mesh.indices[idx_offset].vertex_idx + 1],
                                                       attrib.vertices[3 * shape.mesh.indices[idx_offset].vertex_idx + 2]);

        spt::vec3 raw_normal_l = raw_normal(shape.mesh.indices.begin() + idx_offset, attrib);

        spt::vec3::real_t dot_val = spt::vec3::dot(v0_to_opp, raw_normal_l);
        if (dot_val > 0.0)
            #pragma omp critical
            indices_to_swap.push_back({ shape.mesh.indices.begin() + idx_offset + 1,
                                        shape.mesh.indices.begin() + idx_offset + 2 });
    }

    for (auto& [idx0, idx1] : indices_to_swap)
        std::swap(*idx0, *idx1);
}


void write_file(const tinyobj::attrib_t& attrib, const tinyobj::shape_t& shape,
                std::vector<size_t>& dont_write_faces, std::string_view filename)
{
    std::ofstream file(filename.data());

    for (std::size_t i = 0; i < attrib.vertices.size(); i += 3)
    {
        file << "v "
             << attrib.vertices[i + 0] << ' '
             << attrib.vertices[i + 1] << ' '
             << attrib.vertices[i + 2] << '\n';
    }

    for (std::size_t idx_offset = 0; idx_offset < shape.mesh.indices.size(); idx_offset += 3)
    {
        if (std::find(dont_write_faces.begin(), dont_write_faces.end(), idx_offset) != dont_write_faces.end())
            continue;

        file << "f "
             << shape.mesh.indices[idx_offset + 0].vertex_idx + 1 << ' '
             << shape.mesh.indices[idx_offset + 1].vertex_idx + 1 << ' '
             << shape.mesh.indices[idx_offset + 2].vertex_idx + 1 << '\n';
    }
}


int main()
{
    auto filename = "phset_64_nph_171_phfe.obj";

    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;
    std::string err;
    auto start = std::chrono::steady_clock::now();
    bool ret = tinyobj::load_obj(&attrib, &shapes, &materials, &warn, &err, filename);
    auto stop = std::chrono::steady_clock::now();
    std::chrono::duration<double> elapsed = stop - start;
    std::cout << "Loading time: " << elapsed.count() << " s\n";

    if (!err.empty())
        std::cerr << err << std::endl;

    if (!ret)
        return 1;

    std::vector<size_t> dont_write_faces;
    start = std::chrono::steady_clock::now();
    orientate_mesh(attrib, shapes[0], dont_write_faces);
    stop = std::chrono::steady_clock::now();
    elapsed = stop - start;
    std::cout << "Orientating time: " << elapsed.count() << " s\n";

    write_file(attrib, shapes[0], dont_write_faces, filename);
    std::cout << "done";

    return 0;
}
