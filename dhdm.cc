#include "dhdm.hh"

#include <cassert>
#include <fstream>
#include <iostream>

#include <fmt/core.h>

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

Dhdm::Level::Displacement Dhdm::Level::Iterator::operator *()
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

void Dhdm::Level::Iterator::operator ++()
{
    assert(nrLeft);
    nrLeft--;
}
