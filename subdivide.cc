#include "mesh.hh"
#include "dhdm.hh"

#include <iostream>
#include <memory>

#include <fmt/core.h>

#include <glm/gtx/normal.hpp>

#include <opensubdiv/far/topologyDescriptor.h>
#include <opensubdiv/far/primvarRefiner.h>
#include <opensubdiv/far/stencilTable.h>
#include <opensubdiv/far/stencilTableFactory.h>

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
        if (!dhdms.empty())
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
}
