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

struct Vertex
{
    glm::dvec3 pos;
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
};

double sqr(double x)
{
    return x * x;
}

Mesh readObj(const std::string & path)
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

        std::cerr << fmt::format("rendering tile {}...\n", tile);

        glDrawArrays(GL_TRIANGLES, 0, diff.triangles.size() * 3);

        assert(glGetError() == GL_NO_ERROR);

        std::vector<Pixel> data(width * width);
        glReadBuffer(GL_COLOR_ATTACHMENT0);
        glReadnPixels(0, 0, width, width, GL_RGB, GL_UNSIGNED_SHORT, data.size() * sizeof(Pixel), data.data());

        if (auto err = glGetError())
            throw std::runtime_error("GL error: " + std::to_string(err));

        asyncs.push_back(std::async(std::launch::async, writePNG, tile, std::move(data)));
    }

    std::cerr << fmt::format("waiting...\n");
    for (auto & async : asyncs)
        async.get();
}

int main(int argc, char * * argv)
{
    try {
        if (argc < 4)
            throw std::runtime_error("syntax: dhdm out-prefix base.obj <morph.obj[=weight]...>");

        std::string outPrefix = argv[1];

        auto baseMeshFut = std::async(std::launch::async, &readObj, argv[2]);
        std::vector<std::pair<double, std::future<Mesh>>> morphMeshFuts;
        for (auto i = 3; i < argc; ++i) {
            std::string fn = argv[i];
            double weight = 1.0;
            auto c = fn.find('=');
            if (c != fn.npos) {
                weight = stod(fn.substr(c + 1));
                fn = std::string(fn, 0, c);
            }
            morphMeshFuts.push_back({weight, std::async(std::launch::async, &readObj, fn)});
        }
        auto baseMesh = baseMeshFut.get();
        std::vector<std::pair<double, Mesh>> morphMeshes;
        for (auto & [weight, fut] : morphMeshFuts)
            morphMeshes.push_back({weight, fut.get()});

        auto diff = diffMeshes(1.0, baseMesh, morphMeshes);

        writeDiff(diff, outPrefix);

        return 0;
    } catch (std::exception & e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
