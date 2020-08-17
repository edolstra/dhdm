#pragma once

#include <string>
#include <vector>

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

            Displacement operator *();

            void operator ++();
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

