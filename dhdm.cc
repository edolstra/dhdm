#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <optional>
#include <cassert>
#include <cmath>
#include <future>
#include <set>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <glm/vec3.hpp>
#include <glm/gtx/normal.hpp>
#include <glm/gtx/io.hpp>
#include <fmt/core.h>

#include <png.h>

#include <nlohmann/json.hpp>

#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/stencilTable.h>
#include <opensubdiv/far/stencilTableFactory.h>

#include <boost/iostreams/filtering_stream.hpp>
#include <boost/iostreams/filter/gzip.hpp>

struct Dhdm
{
    std::string contents;

    uint32_t nr_faces = 0;

    struct Level
    {
        uint32_t level;
        uint32_t nrDisplacements;
        std::string_view data;

        struct Displacement
        {
            uint32_t faceIdx;
            uint32_t subfaceIdx;
            unsigned char vertexIdx;
            float x, y, z;
        };

        struct Iterator
        {
            uint32_t level;
            uint32_t nrLeft;
            std::string_view data;
            uint32_t curFaceIdx = 0;
            uint32_t verticesLeft = 0;

            bool operator !=(const Iterator & other) const
            {
                return nrLeft != other.nrLeft;
            }

            Displacement operator *()
            {
                assert(nrLeft);

                if (!verticesLeft) {
                    struct Header
                    {
                        uint32_t faceIdx;
                        uint32_t vertices;
                    };

                    if (data.size() < sizeof(Header))
                        throw std::runtime_error("missing face header");

                    auto header = (Header *) data.data();
                    curFaceIdx = header->faceIdx;
                    verticesLeft = header->vertices;
                    assert(verticesLeft);

                    data = data.substr(sizeof(Header));
                }

                verticesLeft--;

                Displacement displ { .faceIdx = curFaceIdx };
                uint32_t idx;

                if (level < 4) {
                    struct Item
                    {
                        float x;
                        uint8_t b1;
                        uint8_t b2;
                        float y;
                        float z;
                    } __attribute__((packed));

                    static_assert(sizeof(Item) == 14);

                    if (data.size() < sizeof(Item))
                        throw std::runtime_error("missing item");

                    auto item = (Item *) data.data();
                    data = data.substr(sizeof(Item));

                    assert(item->b2 == (level + 1) * 16);

                    displ.x = item->x;
                    displ.y = item->y;
                    displ.z = item->z;
                    idx = ((uint32_t) item->b1) << 8;
                } else {
                    struct Item
                    {
                        float x;
                        uint8_t b1;
                        uint8_t b2;
                        uint8_t b3;
                        uint8_t b4;
                        float y;
                        float z;
                    } __attribute__((packed));

                    static_assert(sizeof(Item) == 16);

                    if (data.size() < sizeof(Item))
                        throw std::runtime_error("missing item");

                    auto item = (Item *) data.data();
                    data = data.substr(sizeof(Item));

                    assert(item->b1 == 0);
                    assert(item->b4 == (level + 1) * 16);

                    displ.x = item->x;
                    displ.y = item->y;
                    displ.z = item->z;
                    idx = ((uint32_t) item->b3) << 8 | ((uint32_t) item->b2);
                }

                displ.subfaceIdx = idx >> (16 - level * 2);
                displ.vertexIdx = (idx >> (14 - level * 2)) & 3;

                return displ;
            }

            void operator ++()
            {
                assert(nrLeft);
                nrLeft--;
            }
        };

        Iterator begin()
        {
            return {.level = level, .nrLeft = nrDisplacements, .data = data};
        }

        Iterator end()
        {
            return {.nrLeft = 0};
        }
    };

    std::vector<Level> levels;

    static constexpr uint32_t MAGIC1 = 0xd0d0d0d0;
    static constexpr uint32_t MAGIC2 = 0x3f800000;

    Dhdm(const std::string & path);
};

Dhdm::Dhdm(const std::string & path)
{
    std::cerr << fmt::format("Parsing '{}'...\n", path);

    std::ifstream fs(path, std::ios::binary);
    if (!fs) throw std::runtime_error("cannot open file");

    contents = std::string(std::istreambuf_iterator<char>(fs), {});
    if (!fs) throw std::runtime_error("cannot read file");

    struct FileHeader
    {
        uint32_t magic1;
        uint32_t nr_levels;
        uint32_t magic2;
        uint32_t nr_levels2;
    };

    if (contents.size() < sizeof(FileHeader))
        throw std::runtime_error("missing file header");

    auto header = (FileHeader *) contents.data();

    if (header->magic1 != MAGIC1 || header->magic2 != MAGIC2)
        throw std::runtime_error("invalid magic");

    if (header->nr_levels != header->nr_levels2)
        throw std::runtime_error("inconsistent number of levels");

    size_t pos = sizeof(FileHeader);

    struct LevelHeader
    {
        uint32_t nr_faces;
        uint32_t level;
        uint32_t nrDisplacements;
        uint32_t size;
    };

    for (uint32_t level = 1; level <= header->nr_levels; ++level) {
        if (contents.size() < pos + sizeof(LevelHeader))
            throw std::runtime_error("missing level header");

        auto lheader = (LevelHeader *) (contents.data() + pos);
        pos += 16;

        std::cerr << fmt::format("  Level {}: {} displacements\n", lheader->level, lheader->nrDisplacements);

        if (level != lheader->level)
            throw std::runtime_error("wrong level header");

        if (nr_faces && nr_faces != lheader->nr_faces)
            throw std::runtime_error("inconsistent number of faces");

        if (contents.size() < pos + sizeof(lheader->size))
            throw std::runtime_error("missing level data");

        levels.push_back(Level {
            .level = level,
            .nrDisplacements = lheader->nrDisplacements,
            .data = std::string_view(contents).substr(pos, lheader->size)
        });

        pos += lheader->size;
    }
}

struct Vertex
{
    glm::dvec3 pos;

    // Interfaces expected by OSD.
    void Clear(void * = nullptr) {
        pos = {0.0, 0.0, 0.0};
    }

    void AddWithWeight(Vertex const & src, float weight) {
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

    void triangulate()
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

    void subdivide(unsigned int level, std::vector<std::pair<double, Dhdm>> dhdms = {});

    static Mesh fromObj(const std::string & path);

    void writeObj(std::ostream & str);

    static Mesh fromDSF(const std::string & geoFile, const std::string & uvFile);
};

void Mesh::subdivide(unsigned int level, std::vector<std::pair<double, Dhdm>> dhdms)
{
    if (level == 0) return;

    std::cerr << fmt::format("Subdividing to level {}...\n", level);

    using namespace OpenSubdiv;

    typedef Far::TopologyDescriptor Descriptor;

    Sdc::SchemeType type = OpenSubdiv::Sdc::SCHEME_CATMARK;

    Sdc::Options options;
    options.SetVtxBoundaryInterpolation(Sdc::Options::VTX_BOUNDARY_EDGE_ONLY);
    options.SetFVarLinearInterpolation(Sdc::Options::FVAR_LINEAR_CORNERS_ONLY);

    std::vector<int> vertsPerFace;
    std::vector<int> vertIndices;
    std::vector<int> uvIndices;
    std::vector<glm::dmat3x3> mats;

    for (auto & face : faces) {
        vertsPerFace.push_back(face.vertices.size());
        for (auto & vert : face.vertices) {
            vertIndices.push_back(vert.vertex);
            uvIndices.push_back(vert.uv);
        }

        assert(face.vertices.size() == 4);

        auto x_axis = glm::normalize(vertices[face.vertices[3].vertex].pos - vertices[face.vertices[0].vertex].pos);
        auto z_axis = glm::normalize(vertices[face.vertices[1].vertex].pos - vertices[face.vertices[0].vertex].pos);
        auto y_axis = -glm::normalize(glm::cross(x_axis, z_axis));
        mats.push_back(glm::dmat3x3(x_axis, y_axis, z_axis));
    }

    Descriptor::FVarChannel uvChannel;
    uvChannel.numValues = uvs.size();
    uvChannel.valueIndices = uvIndices.data();

    Descriptor desc;
    desc.numVertices = vertices.size();
    desc.numFaces = faces.size();
    desc.numVertsPerFace = vertsPerFace.data();
    desc.vertIndicesPerFace = vertIndices.data();
    desc.numFVarChannels = 1;
    desc.fvarChannels = &uvChannel;

    std::unique_ptr<Far::TopologyRefiner> refiner(
        Far::TopologyRefinerFactory<Descriptor>::Create(desc,
            Far::TopologyRefinerFactory<Descriptor>::Options(type, options)));

    {
        Far::TopologyRefiner::UniformOptions refineOptions(level);
        refineOptions.fullTopologyInLastLevel = true;
        refiner->RefineUniform(refineOptions);
    }

    Far::StencilTableFactory::Options stencilOptions;
    stencilOptions.generateIntermediateLevels = true;
    stencilOptions.factorizeIntermediateLevels = false;
    stencilOptions.generateOffsets = true;

    std::unique_ptr<const Far::StencilTable> stencilTable(
        Far::StencilTableFactory::Create(*refiner, stencilOptions));

    std::vector<Vertex> vbuffer(refiner->GetNumVerticesTotal() - vertices.size());

    struct UV
    {
        glm::dvec2 uv;

        // Interfaces expected by OSD.
        void Clear(void * = nullptr) {
            uv = {0.0, 0.0};
        }

        void AddWithWeight(UV const & src, float weight) {
            uv += src.uv * (double) weight;
        }
    };

    std::vector<UV> uvbuffer(refiner->GetNumFVarValuesTotal(0));
    assert(uvs.size() < uvbuffer.size());
    for (size_t i = 0; i < uvs.size(); ++i)
        uvbuffer[i].uv = uvs[i];

    Far::PrimvarRefiner primvarRefiner(*refiner);

    unsigned int start = 0, end = 0;
    UV * srcFVarUV = uvbuffer.data();

    for (unsigned int lvl = 1; lvl <= level; ++lvl) {
        std::cerr << fmt::format("  Applying edits on level {}...\n", lvl);

        auto nverts = refiner->GetLevel(lvl).GetNumVertices();

        auto srcVerts = vertices.data();
        if (lvl > 1) {
             srcVerts = &vbuffer[start];
        }

        start = end;
        end += nverts;

        stencilTable->UpdateValues(srcVerts, vbuffer.data(), start, end);

        const auto & curLevel = refiner->GetLevel(lvl);
        uint32_t nrSubfaces = curLevel.GetNumFaces();

        for (auto & [weight, edits] : dhdms) {
            if (lvl <= edits.levels.size()) {
                for (auto displ : edits.levels[lvl - 1]) {
                    auto subfaceIdx = displ.faceIdx * (1 << (lvl * 2)) + displ.subfaceIdx;
                    assert(subfaceIdx < nrSubfaces);
                    auto fverts = curLevel.GetFaceVertices(subfaceIdx);
                    assert(displ.vertexIdx < fverts.size());
                    auto vertexIdx = fverts[displ.vertexIdx] + start;
                    assert(vertexIdx < vbuffer.size());
                    vbuffer[vertexIdx].pos += weight * mats[displ.faceIdx] * glm::dvec3(displ.x, displ.y, displ.z);
                }
            }
        }

        auto dstFVarUV = srcFVarUV + refiner->GetLevel(lvl - 1).GetNumFVarValues(0);
        primvarRefiner.InterpolateFaceVarying(lvl, srcFVarUV, dstFVarUV, 0);

        srcFVarUV = dstFVarUV;
    }

    vertices.clear();
    faces.clear();
    uvs.clear();

    auto lastLevel = refiner->GetLevel(level);

    auto nverts = lastLevel.GetNumVertices();
    auto nfaces = lastLevel.GetNumFaces();
    auto nuvs   = lastLevel.GetNumFVarValues(0);
    auto firstOfLastVerts = vbuffer.size() - nverts;
    auto firstOfLastUvs = refiner->GetNumFVarValuesTotal(0) - nuvs;

    for (int n = 0; n < nverts; ++n)
        vertices.push_back(vbuffer[firstOfLastVerts + n]);

    for (int n = 0; n < nuvs; ++n)
        uvs.push_back(uvbuffer[firstOfLastUvs + n].uv);

    for (int n = 0; n < nfaces; ++n) {
        std::vector<FaceVertex> vs;
        auto verts = lastLevel.GetFaceVertices(n);
        auto vuvs = lastLevel.GetFaceFVarValues(n, 0);
        assert(verts.size() == vuvs.size());
        for (int i = 0; i < verts.size(); ++i)
            vs.push_back({ .vertex = (VertexId) verts[i], .uv = (UvId) vuvs[i] });
        faces.push_back({ .vertices = std::move(vs) });
    }

    #if 0
    {
        Far::TopologyLevel const & refLastLevel = refiner->GetLevel(level);

        auto nverts = refLastLevel.GetNumVertices();
        auto nfaces = refLastLevel.GetNumFaces();

        auto firstOfLastVerts = vbuffer.size() - nverts;

    }
    #endif
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

    // FIXME: uvs
}

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
};

MeshDiff diffMeshes(
    double scale,
    const Mesh & baseMesh,
    const std::vector<std::pair<double, Mesh>> & morphMeshes)
{
    for (auto & [_, morphMesh] : morphMeshes) {
        assert(baseMesh.vertices.size() == morphMesh.vertices.size());
        assert(baseMesh.faces.size() == morphMesh.faces.size());
    }
    static_assert(sizeof(MeshDiff::GLTriangle) == 3 * 5 * 4);

    /* Compute base mesh vertex normals. */
    std::vector<glm::dvec3> vertexNormals(baseMesh.vertices.size(), glm::dvec3(0.0, 0.0, 0.0));

    for (FaceId faceId = 0; faceId < baseMesh.faces.size(); ++faceId) {
        auto & face = baseMesh.faces[faceId];
        auto normal = glm::triangleNormal(
            baseMesh.vertices[face.vertices[0].vertex].pos,
            baseMesh.vertices[face.vertices[1].vertex].pos,
            baseMesh.vertices[face.vertices[2].vertex].pos);
        vertexNormals[face.vertices[0].vertex] += normal;
        vertexNormals[face.vertices[1].vertex] += normal;
        vertexNormals[face.vertices[2].vertex] += normal;
    }

    for (auto & v : vertexNormals)
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
        return MeshDiff::GLVertex {
            .x = (float) uv.x,
            .y = (float) uv.y,
            .r = (float) (displ2.x * scale + 0.5),
            .g = (float) (displ2.z * scale + 0.5),
            .b = (float) (displ2.y * scale + 0.5),
        };
    };

    for (FaceId faceId = 0; faceId < baseMesh.faces.size(); ++faceId) {
        auto & face = baseMesh.faces[faceId];
        assert(face.vertices.size() == 3);

        /* Compute tangent/bitangent vectors. */
        auto & v0 = baseMesh.vertices[face.vertices[0].vertex].pos;
        auto & v1 = baseMesh.vertices[face.vertices[1].vertex].pos;
        auto & v2 = baseMesh.vertices[face.vertices[2].vertex].pos;

        auto & uv0 = baseMesh.uvs[face.vertices[0].uv];
        auto & uv1 = baseMesh.uvs[face.vertices[1].uv];
        auto & uv2 = baseMesh.uvs[face.vertices[2].uv];

        auto deltaPos1 = v1 - v0;
        auto deltaPos2 = v2 - v0;
        auto deltaUV1 = uv1 - uv0;
        auto deltaUV2 = uv2 - uv0;

        auto r = 1.0 / (deltaUV1.x * deltaUV2.y - deltaUV1.y * deltaUV2.x);
        auto tangent = (deltaPos1 * deltaUV2.y - deltaPos2 * deltaUV1.y) * r;
        auto bitangent = (deltaPos2 * deltaUV1.x - deltaPos1 * deltaUV2.x) * r;

        /* Compute the object space displacements. */
        glm::dvec3 d0, d1, d2;
        for (auto & [weight, morphMesh] : morphMeshes) {
            auto & face2 = morphMesh.faces[faceId];
            d0 += weight * (morphMesh.vertices[face2.vertices[0].vertex].pos - baseMesh.vertices[face.vertices[0].vertex].pos);
            d1 += weight * (morphMesh.vertices[face2.vertices[1].vertex].pos - baseMesh.vertices[face.vertices[1].vertex].pos);
            d2 += weight * (morphMesh.vertices[face2.vertices[2].vertex].pos - baseMesh.vertices[face.vertices[2].vertex].pos);
        }

        diff.triangles.push_back(
            MeshDiff::GLTriangle {
                .v0 = toVertex(tangent, bitangent, vertexNormals[face.vertices[0].vertex], d0, baseMesh.uvs[face.vertices[0].uv]),
                .v1 = toVertex(tangent, bitangent, vertexNormals[face.vertices[1].vertex], d1, baseMesh.uvs[face.vertices[1].uv]),
                .v2 = toVertex(tangent, bitangent, vertexNormals[face.vertices[2].vertex], d2, baseMesh.uvs[face.vertices[2].uv]),
            });
    }

    std::cerr << "maximum displacement: " << maxDispl << "\n";
    std::cerr << "maximum pixel component value: " << (maxDispl * scale + 0.5) * 255.0 << "\n";

    return diff;
}

const char * vertexSource = R"glsl(
#version 150 core

uniform float tile;

in vec2 position;
in vec3 color;

out vec3 Color;

void main()
{
    Color = color;
    gl_Position = vec4(2.0 * (position.x - 0.5 - tile), 2.0 * (position.y - 0.5), 0.0, 1.0);
}
)glsl";

const char* fragmentSource = R"glsl(
#version 150 core

in  vec3 Color;
out vec4 outColor;

void main()
{
    outColor = vec4(Color, 1.0);
}
)glsl";

void writeDiff(const MeshDiff & diff, const std::string & outPrefix)
{
    if (!glfwInit())
        throw std::runtime_error("failed to initialize GLFW");

    glfwWindowHint(GLFW_SAMPLES, 4); // 4x antialiasing
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3); // We want OpenGL 3.3
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // To make MacOS happy; should not be needed
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE); // We don't want the old OpenGL
    glfwWindowHint(GLFW_VISIBLE, GLFW_FALSE);

    uint16_t width = 4096;

    // Open a window and create its OpenGL context
    GLFWwindow * window; // (In the accompanying source code, this variable is global for simplicity)
    window = glfwCreateWindow(width, width, "Tutorial 01", nullptr, nullptr);
    if (!window)
        throw std::runtime_error("failed to open GLFW window");
    glfwMakeContextCurrent(window); // Initialize GLEW
    glewExperimental = true; // Needed in core profile
    if (glewInit() != GLEW_OK)
        throw std::runtime_error("failed to initialize GLEW");

    GLuint vao;
    glGenVertexArrays(1, &vao);
    glBindVertexArray(vao);

    GLuint vertexBuffer;
    glGenBuffers(1, &vertexBuffer);
    glBindBuffer(GL_ARRAY_BUFFER, vertexBuffer);
    glBufferData(GL_ARRAY_BUFFER,
        diff.triangles.size() * sizeof(MeshDiff::GLTriangle),
        diff.triangles.data(), GL_STATIC_DRAW);


    GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
    glShaderSource(vertexShader, 1, &vertexSource, nullptr);
    glCompileShader(vertexShader);
    GLint status;
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &status);
    assert(status == GL_TRUE);

    GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    glShaderSource(fragmentShader, 1, &fragmentSource, nullptr);
    glCompileShader(fragmentShader);
    glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &status);
    assert(status == GL_TRUE);

    GLuint shaderProgram = glCreateProgram();
    glAttachShader(shaderProgram, vertexShader);
    glAttachShader(shaderProgram, fragmentShader);
    glBindFragDataLocation(shaderProgram, 0, "outColor");

    glLinkProgram(shaderProgram);
    glUseProgram(shaderProgram);

    GLint posAttrib = glGetAttribLocation(shaderProgram, "position");
    GLint colAttrib = glGetAttribLocation(shaderProgram, "color");
    GLint tileAttrib = glGetUniformLocation(shaderProgram, "tile");

    glEnableVertexAttribArray(posAttrib);
    glVertexAttribPointer(posAttrib, 2, GL_FLOAT, GL_FALSE, sizeof(MeshDiff::GLVertex), nullptr);

    glEnableVertexAttribArray(colAttrib);
    glVertexAttribPointer(colAttrib, 3, GL_FLOAT, GL_FALSE, sizeof(MeshDiff::GLVertex), (void*) (2 * sizeof(float)));


    GLuint frameBuffer;
    glGenFramebuffers(1, &frameBuffer);

    glBindFramebuffer(GL_FRAMEBUFFER, frameBuffer);

    GLuint texColorBuffer;
    glGenTextures(1, &texColorBuffer);
    glBindTexture(GL_TEXTURE_2D, texColorBuffer);

    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB16F, width, width, 0, GL_RGBA, GL_UNSIGNED_SHORT, nullptr);
    assert(glGetError() == GL_NO_ERROR);

    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, texColorBuffer, 0);


    struct Pixel
    {
        uint16_t r, g, b;
    };

    static_assert(sizeof(Pixel) == 6);

    auto writePNG = [&](unsigned int tile, std::vector<Pixel> data)
    {
        /* Blur uninitialized pixels by averaging adjacent initialized
           pixels. */
        for (int x = 0; x < width; x++)
            for (int y = 0; y < width; y++) {
                auto & pixel = data[y * width + x];
                if (pixel.r == 0 && pixel.g == 0 && pixel.b == 0) {
                    uint32_t r = 0, g = 0, b = 0;
                    unsigned int c = 0;
                    auto read = [&](int x, int y)
                    {
                        if (x >= 0 && x < width && y >= 0 && y < width) {
                            auto & pixel2 = data[y * width + x];
                            if (pixel2.r != 0 || pixel2.g != 0 || pixel2.b != 0) {
                                r += pixel2.r;
                                g += pixel2.g;
                                b += pixel2.b;
                                c++;
                            }
                        }
                    };
                    read(x - 1, y - 1);
                    read(x,     y - 1);
                    read(x + 1, y - 1);
                    read(x - 1, y);
                    read(x + 1, y);
                    read(x - 1, y + 1);
                    read(x,     y + 1);
                    read(x + 1, y + 1);
                    if (c) {
                        //std::cerr << fmt::format("FIX {} {} {}", x, y, c);
                        pixel.r = r / c;
                        pixel.g = g / c;
                        pixel.b = b / c;
                    }
                }
            }

        for (int x = 0; x < width; x++)
            for (int y = 0; y < width; y++) {
                auto & pixel = data[y * width + x];
                if (pixel.r == 0 && pixel.g == 0 && pixel.b == 0) {
                    pixel.r = 32768;
                    pixel.g = 32768;
                    pixel.b = 32768;
                }
            }

        /* Write the PNG. */
        auto outFile = fmt::format("{}-{}.png", outPrefix, tile);
        std::cerr << fmt::format("writing '{}'...\n", outFile);

        auto fp = fopen(outFile.c_str(), "wb");
        assert(fp);

        auto png_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, nullptr, nullptr, nullptr);
        assert(png_ptr);

        auto info_ptr = png_create_info_struct(png_ptr);
        assert(info_ptr);

        png_init_io(png_ptr, fp);

        png_set_IHDR(
            png_ptr, info_ptr, width, width,
            16, PNG_COLOR_TYPE_RGB, PNG_INTERLACE_NONE,
            PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

        png_write_info(png_ptr, info_ptr);

        png_byte * row_pointers[width];
        for (size_t i = 0; i < width; i++)
            row_pointers[i] = ((png_byte *) data.data()) + (width - i - 1) * width * 6;

        png_set_swap(png_ptr);

        png_write_image(png_ptr, row_pointers);

        png_write_end(png_ptr, nullptr);

        png_destroy_read_struct(&png_ptr, &info_ptr, nullptr);
    };


    std::vector<std::future<void>> asyncs;

    for (auto & tile : diff.tiles) {

        glUniform1f(tileAttrib, tile);

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        std::cerr << fmt::format("Rendering tile {}...\n", tile);

        glDrawArrays(GL_TRIANGLES, 0, diff.triangles.size() * 3);

        assert(glGetError() == GL_NO_ERROR);

        std::vector<Pixel> data(width * width);
        glReadBuffer(GL_COLOR_ATTACHMENT0);
        glReadnPixels(0, 0, width, width, GL_RGB, GL_UNSIGNED_SHORT, data.size() * sizeof(Pixel), data.data());

        if (auto err = glGetError())
            throw std::runtime_error("GL error: " + std::to_string(err));

        asyncs.push_back(std::async(std::launch::async, writePNG, tile, std::move(data)));
    }

    std::cerr << fmt::format("Waiting...\n");
    for (auto & async : asyncs)
        async.get();
}

void applyMorphs(Mesh & mesh, size_t nrPaths, char * paths[])
{
    std::vector<std::pair<double, Dhdm>> dhdms;

    for (size_t i = 0; i < nrPaths; ++i) {
        std::string fn = paths[i];
        double weight = 1.0;
        auto c = fn.find('=');
        if (c != fn.npos) {
            weight = stod(fn.substr(c + 1));
            fn = std::string(fn, 0, c);
        }

        if (fn.size() >= 4 && std::string(fn, fn.size() - 4) == ".dsf") {
            std::cerr << fmt::format("Applying morph '{}'...\n", fn);
            auto morph = readJSON(fn)["modifier_library"][0]["morph"];
            int32_t vertexCount = morph["vertex_count"];
            if (vertexCount != -1 && (size_t) vertexCount != mesh.vertices.size())
                throw std::runtime_error(
                    fmt::format("Morph vertex count {} differs from base mesh vertex count {}.",
                        vertexCount,
                        mesh.vertices.size()));
            for (auto v : morph["deltas"]["values"])
                mesh.vertices[v[0]].pos += weight * glm::dvec3(v[1], -(double) v[3], v[2]);
            if (morph.find("hd_url") != morph.end()) {
                dhdms.push_back({weight, Dhdm(fn.substr(0, fn.size() - 4) + ".dhdm")});
                //dhdms.push_back({weight, Dhdm(percentDecode((std::string) morph["hd_url"]))});
            }
        } else {
            dhdms.push_back({weight, Dhdm(fn)});
        }
    }

    mesh.subdivide(4, dhdms);
}

int main(int argc, char * * argv)
{
    try {

        if (argc < 2)
            throw std::runtime_error(
                "Syntax: dhdm morphs-to-displacement <prefix> <base-mesh>.dsf <uv-map>.dsf [<morph>.[dsf|dhdm][=weight]]...\n"
                "           | morphs-to-obj <base-mesh>.dsf <uv-map>.dsf [<morph>.[dsf|dhdm][=weight]]... > dest.obj\n"
                "           | diff-objs <prefix> <base-mesh>.obj [<mesh>.obj[=weight]]...\n");

        std::string verb = argv[1];

        if (verb == "morphs-to-displacement") {
            if (argc < 5)
                throw std::runtime_error("Missing arguments.");

            std::string prefix = argv[2];
            std::string baseMeshPath = argv[3];
            std::string uvMapPath = argv[4];

            auto baseMesh = Mesh::fromDSF(baseMeshPath, uvMapPath);
            auto finalMesh = baseMesh;
            applyMorphs(finalMesh, argc - 5, argv + 5);

            baseMesh.subdivide(4, {});

            baseMesh.triangulate();
            finalMesh.triangulate();

            auto diff = diffMeshes(1.0, baseMesh, {{1.0, std::move(finalMesh)}});

            writeDiff(diff, prefix);
        }

        else if (verb == "morphs-to-obj") {
            if (argc < 4)
                throw std::runtime_error("Missing arguments.");

            std::string baseMeshPath = argv[2];
            std::string uvMapPath = argv[3];

            auto baseMesh = Mesh::fromDSF(baseMeshPath, uvMapPath);
            applyMorphs(baseMesh, argc - 4, argv + 4);
            baseMesh.writeObj(std::cout);
        }

        else if (verb == "diff-objs") {

            std::string prefix = argv[2];
            std::string baseMeshPath = argv[3];

            auto baseMeshFut = std::async(std::launch::async, &Mesh::fromObj, baseMeshPath);
            std::vector<std::pair<double, std::future<Mesh>>> morphMeshFuts;
            for (auto i = 4; i < argc; ++i) {
                std::string fn = argv[i];
                double weight = 1.0;
                auto c = fn.find('=');
                if (c != fn.npos) {
                    weight = stod(fn.substr(c + 1));
                    fn = std::string(fn, 0, c);
                }
                morphMeshFuts.push_back({weight, std::async(std::launch::async, &Mesh::fromObj, fn)});
            }
            auto baseMesh = baseMeshFut.get();
            std::vector<std::pair<double, Mesh>> morphMeshes;
            for (auto & [weight, fut] : morphMeshFuts)
                morphMeshes.push_back({weight, fut.get()});

            auto diff = diffMeshes(1.0, baseMesh, morphMeshes);

            writeDiff(diff, prefix);
        }

        else
            throw std::runtime_error(fmt::format("Unknown command '{}'.", verb));

        return 0;
    } catch (std::exception & e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
