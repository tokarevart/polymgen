#include <iostream>
#include <algorithm>
#include <string>
#include <ctime>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include "vec.h"


using indices = std::vector<tinyobj::index_t>;


#define TIME_TO_MS(time_t_var) \
    (time_t_var * 1000 / CLOCKS_PER_SEC)


bool faceContainsVert( indices::const_iterator face, indices::const_iterator vert )
{
    if (   (face + 0)->vertex_idx == vert->vertex_idx
        || (face + 1)->vertex_idx == vert->vertex_idx
        || (face + 2)->vertex_idx == vert->vertex_idx)
        return true;

    return false;
}


bool findOneVertNot( int& out, indices::const_iterator inFace, indices::const_iterator notFace )
{
    int n_not_conts = 0;
    bool not_conts[3];
    if (not_conts[0] = !faceContainsVert(notFace, inFace + 0)) n_not_conts++;
    if (not_conts[1] = !faceContainsVert(notFace, inFace + 1)) n_not_conts++;
    if (not_conts[2] = !faceContainsVert(notFace, inFace + 2)) n_not_conts++;
    if (n_not_conts != 1) return false;

    if (not_conts[0])
    {
        out = (inFace + 0)->vertex_idx;
        return true;
    }
    if (not_conts[1])
    {
        out = (inFace + 1)->vertex_idx;
        return true;
    }
    if (not_conts[2])
    {
        out = (inFace + 2)->vertex_idx;
        return true;
    }

    throw std::logic_error("Function findOneVertNot didn't find vertex.");
}


int findOppVert( indices::const_iterator face, const tinyobj::shape_t& shape )
{
    std::vector<int> opp_idces;

    for (size_t i = 0; i < shape.mesh.indices.size(); i += 3)
    {
        int opp_idx;
        if (findOneVertNot(opp_idx, shape.mesh.indices.begin() + i, face))
        {
            if (std::count(opp_idces.begin(), opp_idces.end(), opp_idx) == 2)
                return opp_idx;

            opp_idces.push_back(opp_idx);
        }
    }

    return -1;
}


vec3 getVertPos( int idx, const tinyobj::attrib_t& attrib )
{
    return vec3(attrib.vertices[3 * idx + 0],
                attrib.vertices[3 * idx + 1],
                attrib.vertices[3 * idx + 2]);
}


vec3 computeRawNormal( indices::const_iterator face, const tinyobj::attrib_t& attrib )
{
    vec3 v0(attrib.vertices[3 * (face + 0)->vertex_idx + 0],
            attrib.vertices[3 * (face + 0)->vertex_idx + 1],
            attrib.vertices[3 * (face + 0)->vertex_idx + 2]);
    vec3 v1(attrib.vertices[3 * (face + 1)->vertex_idx + 0],
            attrib.vertices[3 * (face + 1)->vertex_idx + 1],
            attrib.vertices[3 * (face + 1)->vertex_idx + 2]);
    vec3 v2(attrib.vertices[3 * (face + 2)->vertex_idx + 0],
            attrib.vertices[3 * (face + 2)->vertex_idx + 1],
            attrib.vertices[3 * (face + 2)->vertex_idx + 2]);

    return vec3::cross(v1 - v0, v2 - v0);
}


void orientateMesh( const tinyobj::attrib_t& attrib, tinyobj::shape_t& shape, std::vector<size_t>& dontWriteFaces )
{
    std::vector<std::pair<indices::iterator, indices::iterator>> indices_to_swap;

    size_t n_indices = shape.mesh.indices.size();
    #pragma omp parallel for
    for (size_t idx_offset = 0; idx_offset < n_indices; idx_offset += 3)
    {
        int opp_vert_idx = findOppVert(shape.mesh.indices.begin() + idx_offset, shape);
        if (opp_vert_idx == -1)
            #pragma omp critical
            dontWriteFaces.push_back(idx_offset);
        vec3 opp_vert_pos = getVertPos(opp_vert_idx, attrib);

        vec3 v0_to_opp = opp_vert_pos - vec3(attrib.vertices[3 * shape.mesh.indices[idx_offset].vertex_idx + 0],
                                             attrib.vertices[3 * shape.mesh.indices[idx_offset].vertex_idx + 1],
                                             attrib.vertices[3 * shape.mesh.indices[idx_offset].vertex_idx + 2]);

        vec3 raw_normal = computeRawNormal(shape.mesh.indices.begin() + idx_offset, attrib);

        float dot_val = vec3::dot(v0_to_opp, raw_normal);
        if (dot_val > 0.0f)
            #pragma omp critical
            indices_to_swap.push_back({ shape.mesh.indices.begin() + idx_offset + 1,
                                        shape.mesh.indices.begin() + idx_offset + 2 });
    }

    for (auto& idces : indices_to_swap)
        std::swap(*idces.first, *idces.second);
}


void writeFile( const tinyobj::attrib_t& attrib, const tinyobj::shape_t& shape, std::vector<size_t>& dontWriteFaces, std::string_view filename )
{
    std::ofstream file(filename.data());

    for (size_t i = 0; i < attrib.vertices.size(); i += 3)
    {
        file << "v "
             << attrib.vertices[i + 0] << ' '
             << attrib.vertices[i + 1] << ' '
             << attrib.vertices[i + 2] << '\n';
    }

    for (size_t idx_offset = 0; idx_offset < shape.mesh.indices.size(); idx_offset += 3)
    {
        if (std::find(dontWriteFaces.begin(), dontWriteFaces.end(), idx_offset) != dontWriteFaces.end())
            continue;

        file << "f "
             << shape.mesh.indices[idx_offset + 0].vertex_idx + 1 << ' '
             << shape.mesh.indices[idx_offset + 1].vertex_idx + 1 << ' '
             << shape.mesh.indices[idx_offset + 2].vertex_idx + 1 << '\n';
    }
}


int main()
{
    tinyobj::attrib_t attrib;
    std::vector<tinyobj::shape_t> shapes;
    std::vector<tinyobj::material_t> materials;

    std::string warn;
    std::string err;
    time_t start = clock();
    bool ret = tinyobj::load_obj(&attrib, &shapes, &materials, &warn, &err, "phset.obj");
    time_t elapsed = clock() - start;
    std::cout << "Loading time: " << TIME_TO_MS(elapsed) << " ms\n";

    if (!err.empty())
        std::cerr << err << std::endl;

    if (!ret)
        return 1;

    std::vector<size_t> dont_write_faces;
    start = clock();
    orientateMesh(attrib, shapes[0], dont_write_faces);
    elapsed = clock() - start;
    std::cout << "Orientating time: " << TIME_TO_MS(elapsed) << " ms\n";

    writeFile(attrib, shapes[0], dont_write_faces, "phset.obj");
    std::cout << "done";

    return 0;
}
