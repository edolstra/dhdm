#include "mesh.hh"

#include <optional>
#include <fstream>
#include <iostream>

#include <glm/gtx/io.hpp>

#include <fmt/core.h>

void Mesh::triangulate()
{
    auto nrFaces = faces.size();
    faces.reserve(nrFaces * 2);
    for (size_t i = 0; i < nrFaces; ++i) {
        auto & face = faces[i];
        if (face.vertices.size() == 3) continue;
        assert(face.vertices.size() == 4);
        faces.push_back(Face { .vertices = { face.vertices[2], face.vertices[3], face.vertices[0] } });
        face.vertices.pop_back();
    }
}

Mesh Mesh::fromObj(const std::string & path)
{
    auto fs = std::fstream(path, std::fstream::in);
    if (!fs) throw std::runtime_error("cannot open file");

    Mesh mesh;

    std::optional<std::string> objectName;

    std::string line;
    while (std::getline(fs, line)) {
        if (std::string_view(line).substr(0, 2) == "o ") {
            if (objectName)
                throw std::runtime_error("obj file contains multiple meshes");
            objectName = std::string(line, 2);
            std::cerr << "mesh name: " << *objectName << "\n";
        } else if (std::string_view(line).substr(0, 2) == "v ") {
            double x, y, z;
            if (sscanf(line.c_str() + 2, "%lf %lf %lf", &x, &y, &z) != 3)
                throw std::runtime_error("invalid vertex: " + line);
            mesh.vertices.push_back({glm::dvec3(x, -z, y)});
        } else if (std::string_view(line).substr(0, 3) == "vt ") {
            double u, v;
            if (sscanf(line.c_str() + 3, "%lf %lf", &u, &v) != 2)
                throw std::runtime_error("invalid texture coordinate: " + line);
            mesh.uvs.push_back({ u, v });
        } else if (std::string_view(line).substr(0, 2) == "f ") {
            uint32_t v1, v2, v3, v4;
            uint32_t uv1, uv2, uv3, uv4;
            if (sscanf(line.c_str() + 2, "%d/%d %d/%d %d/%d %d/%d", &v1, &uv1, &v2, &uv2, &v3, &uv3, &v4, &uv4) != 8)
                throw std::runtime_error("invalid face: " + line);
            assert(v1 >= 1 && v1 <= mesh.vertices.size());
            assert(v2 >= 1 && v2 <= mesh.vertices.size());
            assert(v3 >= 1 && v3 <= mesh.vertices.size());
            assert(v4 >= 1 && v4 <= mesh.vertices.size());
            assert(uv1 >= 1 && uv1 <= mesh.uvs.size());
            assert(uv2 >= 1 && uv2 <= mesh.uvs.size());
            assert(uv3 >= 1 && uv3 <= mesh.uvs.size());
            assert(uv4 >= 1 && uv4 <= mesh.uvs.size());
            Face face;
            face.vertices.push_back({ .vertex = VertexId(v1 - 1), .uv = UvId(uv1 - 1) });
            face.vertices.push_back({ .vertex = VertexId(v2 - 1), .uv = UvId(uv2 - 1) });
            face.vertices.push_back({ .vertex = VertexId(v3 - 1), .uv = UvId(uv3 - 1) });
            face.vertices.push_back({ .vertex = VertexId(v4 - 1), .uv = UvId(uv4 - 1) });
            mesh.faces.push_back(std::move(face));
        }
    }

    std::cerr << "number of vertices: " << mesh.vertices.size() << "\n";
    std::cerr << "number of faces: " << mesh.faces.size() << "\n";

    mesh.triangulate();

    std::cerr << "number of tris: " << mesh.faces.size() << "\n";

    return mesh;
}

void Mesh::writeObj(std::ostream & str)
{
    std::cerr << "Writing obj...\n";

    for (auto & vert : vertices)
        str << fmt::format("v {} {} {}\n", (float) vert.pos.x, (float) vert.pos.y, (float) vert.pos.z);

    for (auto & uv : uvs)
        str << fmt::format("vt {} {}\n", (float) uv.x, (float) uv.y);

    for (auto & face : faces) {
        str << "f ";
        for (auto & vert : face.vertices)
            str << (vert.vertex + 1) << "/" << (vert.uv + 1) << " ";
        str << "\n";
    }
}
