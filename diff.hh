#pragma once

#include <vector>
#include <set>

struct Mesh;

struct MeshDiff
{
    struct GLVertex
    {
        float x, y;
        float r, g, b;
    };

    struct GLTriangle
    {
        GLVertex v0, v1, v2;
    };

    std::vector<GLTriangle> triangles;

    std::set<unsigned int> tiles;

    static MeshDiff create(
        double scale,
        const Mesh & baseMesh,
        const std::vector<std::pair<double, Mesh>> & morphMeshes);

    void writeToPNG(const std::string & outPrefix);
};
