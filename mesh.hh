#pragma once

#include <optional>
#include <string>
#include <vector>

#include <glm/vec2.hpp>
#include <glm/vec3.hpp>

struct Dhdm;

struct Vertex
{
    glm::dvec3 pos;

    // Interfaces expected by OSD.
    void Clear(void * = nullptr)
    {
        pos = {0.0, 0.0, 0.0};
    }

    void AddWithWeight(Vertex const & src, float weight)
    {
        pos += src.pos * (double) weight;
    }
};

typedef uint32_t VertexId;
typedef uint32_t UvId;
typedef uint32_t FaceId;

struct FaceVertex
{
    VertexId vertex;
    UvId uv;
};

struct Face
{
    std::vector<FaceVertex> vertices;
};

struct Mesh
{
    std::vector<Vertex> vertices;
    std::vector<glm::dvec2> uvs;
    std::vector<Face> faces;

    void triangulate();

    void subdivide(unsigned int level, std::vector<std::pair<double, Dhdm>> dhdms = {});

    static Mesh fromObj(const std::string & path);

    void writeObj(std::ostream & str);

    static Mesh fromDSF(const std::string & geoFile, const std::string & uvFile);

    std::optional<std::string> applyMorph(double weight, const std::string & dsfFile);
};
