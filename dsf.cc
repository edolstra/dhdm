#include "mesh.hh"

#include <iostream>
#include <fstream>

#include <nlohmann/json.hpp>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#include <fmt/core.h>

nlohmann::json readJSON(const std::string & path)
{
    auto fs = std::fstream(path, std::fstream::in);
    if (!fs) throw std::runtime_error("cannot open file");

    unsigned char magic[2];
    fs.read((char *) magic, sizeof(magic));
    fs.seekg(0, fs.beg);

    if (magic[0] == 0x1f && magic[1] == 0x8b) {
        boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
        in.push(boost::iostreams::gzip_decompressor());
        in.push(fs);
        std::istream str(&in);
        nlohmann::json json;
        str >> json;
        return json;
    } else {
        nlohmann::json json;
        fs >> json;
        return json;
    }
}

template <class T>
std::optional<typename T::mapped_type> get(const T & map, const typename T::key_type & key)
{
    auto i = map.find(key);
    if (i == map.end()) return {};
    return std::optional<typename T::mapped_type>(i->second);
}

Mesh Mesh::fromDSF(const std::string & geoFile, const std::string & uvFile)
{
    Mesh mesh;

    std::map<std::pair<FaceId, VertexId>, UvId> overrides;

    {
        auto uvMap = readJSON(uvFile)["uv_set_library"][0];
        for (auto & uv : uvMap["uvs"]["values"]) {
            assert(uv.size() == 2);
            mesh.uvs.push_back({glm::dvec2(uv[0], uv[1])});
        }
        for (auto & p : uvMap["polygon_vertex_indices"]) {
            assert(p.size() == 3);
            overrides.insert_or_assign({p[0], p[1]}, p[2]);
        }
    }

    auto geometry = readJSON(geoFile)["geometry_library"][0];

    for (auto & vertex : geometry["vertices"]["values"]) {
        assert(vertex.size() == 3);
        mesh.vertices.push_back({glm::dvec3(vertex[0], -(double) vertex[2], vertex[1])});
    }

    for (auto & poly : geometry["polylist"]["values"]) {
        assert(poly.size() >= 5);
        std::vector<FaceVertex> vertices;
        FaceId faceIdx = mesh.faces.size();
        for (size_t i = 2; i < poly.size(); ++i) {
            VertexId vertexIdx = poly[i];
            assert(vertexIdx < mesh.vertices.size());
            auto uvIdx = get(overrides, {faceIdx, vertexIdx}).value_or(vertexIdx);
            vertices.push_back({ .vertex = vertexIdx, .uv = uvIdx });
        }
        mesh.faces.push_back({ .vertices = std::move(vertices) });
    }

    std::cerr << fmt::format("Read {}: {} vertices, {} faces, {} UVs\n",
        geoFile, mesh.vertices.size(), mesh.faces.size(), mesh.uvs.size());

    return mesh;
}

std::optional<std::string> Mesh::applyMorph(double weight, const std::string & dsfFile)
{
    std::cerr << fmt::format("Applying morph '{}'...\n", dsfFile);
    auto morph = readJSON(dsfFile)["modifier_library"][0]["morph"];
    int32_t vertexCount = morph["vertex_count"];
    if (vertexCount != -1 && (size_t) vertexCount != vertices.size())
        throw std::runtime_error(
            fmt::format("Morph vertex count {} differs from base mesh vertex count {}.",
                vertexCount,
                vertices.size()));
    for (auto v : morph["deltas"]["values"])
        vertices[v[0]].pos += weight * glm::dvec3(v[1], -(double) v[3], v[2]);
    if (morph.find("hd_url") != morph.end())
        return dsfFile.substr(0, dsfFile.size() - 4) + ".dhdm";
    else
        return {};
}
