#include "diff.hh"
#include "mesh.hh"

#include <iostream>

#include <glm/gtx/normal.hpp>

MeshDiff MeshDiff::create(
    double scale,
    const Mesh & baseMesh,
    const std::vector<std::pair<double, Mesh>> & morphMeshes)
{
    for (auto & [_, morphMesh] : morphMeshes) {
        assert(baseMesh.vertices.size() == morphMesh.vertices.size());
        assert(baseMesh.faces.size() == morphMesh.faces.size());
    }
    static_assert(sizeof(GLTriangle) == 3 * 5 * 4);

    /* Compute base mesh vertex normals, tangents and bitangents. */
    std::vector<glm::dvec3> vertexNormals(baseMesh.vertices.size(), glm::dvec3(0.0, 0.0, 0.0));
    std::vector<glm::dvec3> vertexTangents(baseMesh.vertices.size(), glm::dvec3(0.0, 0.0, 0.0));
    std::vector<glm::dvec3> vertexBitangents(baseMesh.vertices.size(), glm::dvec3(0.0, 0.0, 0.0));

    for (FaceId faceId = 0; faceId < baseMesh.faces.size(); ++faceId) {
        auto & face = baseMesh.faces[faceId];
        assert(face.vertices.size() == 3);

        auto vi0 = face.vertices[0].vertex;
        auto vi1 = face.vertices[1].vertex;
        auto vi2 = face.vertices[2].vertex;

        auto & v0 = baseMesh.vertices[vi0].pos;
        auto & v1 = baseMesh.vertices[vi1].pos;
        auto & v2 = baseMesh.vertices[vi2].pos;

        auto normal = glm::triangleNormal(v0, v1, v2);
        vertexNormals[vi0] += normal;
        vertexNormals[vi1] += normal;
        vertexNormals[vi2] += normal;

        auto & uv0 = baseMesh.uvs[face.vertices[0].uv];
        auto & uv1 = baseMesh.uvs[face.vertices[1].uv];
        auto & uv2 = baseMesh.uvs[face.vertices[2].uv];

        auto deltaPos1 = v1 - v0;
        auto deltaPos2 = v2 - v0;
        auto deltaUV1 = uv1 - uv0;
        auto deltaUV2 = uv2 - uv0;

        auto r = 1.0 / (deltaUV1.x * deltaUV2.y - deltaUV1.y * deltaUV2.x);

        auto tangent = (deltaPos1 * deltaUV2.y - deltaPos2 * deltaUV1.y) * r;
        vertexTangents[vi0] += tangent;
        vertexTangents[vi1] += tangent;
        vertexTangents[vi2] += tangent;

        auto bitangent = (deltaPos2 * deltaUV1.x - deltaPos1 * deltaUV2.x) * r;
        vertexBitangents[vi0] += bitangent;
        vertexBitangents[vi1] += bitangent;
        vertexBitangents[vi2] += bitangent;
    }

    for (auto & v : vertexNormals)
        v = glm::normalize(v);

    for (auto & v : vertexTangents)
        v = glm::normalize(v);

    for (auto & v : vertexBitangents)
        v = glm::normalize(v);

    MeshDiff diff;

    double maxDispl = 0.0;

    auto toVertex = [&](glm::dvec3 tangent, glm::dvec3 bitangent, const glm::dvec3 & normal, const glm::dvec3 & displ, const glm::dvec2 & uv)
    {
        // Gram-Schmidt orthogonalization.
        tangent = glm::normalize(tangent - normal * glm::dot(normal, tangent));
        bitangent = glm::normalize(bitangent - normal * glm::dot(normal, bitangent) - tangent * glm::dot(tangent, bitangent));

        if (glm::dot(glm::cross(normal, tangent), bitangent) < 0.0)
            tangent = tangent * -1.0;

        // Convert displacement from object to normal space.
        auto mat = glm::dmat3x3(
            tangent.x, bitangent.x, normal.x,
            tangent.y, bitangent.y, normal.y,
            tangent.z, bitangent.z, normal.z);

        glm::dvec3 displ2 = mat * displ;

        auto max = std::max(abs(displ2.x), std::max(abs(displ2.y), abs(displ2.z)));
        maxDispl = std::max(maxDispl, max);

        if (max > 1e-6) {
            unsigned int tile = uv.x / 1.0;
            diff.tiles.insert(tile);
        }

        // Note: Blender expects the displacement along the normal
        // axis to be in the green channel, the tangent axis in the
        // red channel, and the bitangent axis in the blue channel.
        return GLVertex {
            .x = (float) uv.x,
            .y = (float) uv.y,
            .r = (float) (displ2.x * scale + 0.5),
            .g = (float) (displ2.z * scale + 0.5),
            .b = (float) (displ2.y * scale + 0.5),
        };
    };

    for (FaceId faceId = 0; faceId < baseMesh.faces.size(); ++faceId) {
        auto & face = baseMesh.faces[faceId];

        auto vi0 = face.vertices[0].vertex;
        auto vi1 = face.vertices[1].vertex;
        auto vi2 = face.vertices[2].vertex;

        /* Compute the object space displacements. */
        glm::dvec3 d0, d1, d2;
        for (auto & [weight, morphMesh] : morphMeshes) {
            auto & face2 = morphMesh.faces[faceId];
            d0 += weight * (morphMesh.vertices[face2.vertices[0].vertex].pos - baseMesh.vertices[vi0].pos);
            d1 += weight * (morphMesh.vertices[face2.vertices[1].vertex].pos - baseMesh.vertices[vi1].pos);
            d2 += weight * (morphMesh.vertices[face2.vertices[2].vertex].pos - baseMesh.vertices[vi2].pos);
        }

        diff.triangles.push_back(
            GLTriangle {
                .v0 = toVertex(vertexTangents[vi0], vertexBitangents[vi0], vertexNormals[vi0], d0, baseMesh.uvs[face.vertices[0].uv]),
                .v1 = toVertex(vertexTangents[vi1], vertexBitangents[vi1], vertexNormals[vi1], d1, baseMesh.uvs[face.vertices[1].uv]),
                .v2 = toVertex(vertexTangents[vi2], vertexBitangents[vi2], vertexNormals[vi2], d2, baseMesh.uvs[face.vertices[2].uv]),
            });
    }

    std::cerr << "Maximum displacement: " << maxDispl << "\n";
    std::cerr << "Maximum pixel component value: " << (maxDispl * scale + 0.5) * 255.0 << "\n";

    return diff;
}
