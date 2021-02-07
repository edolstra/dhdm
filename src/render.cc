#include "diff.hh"

#include <cassert>
#include <future>
#include <iostream>

#include <GL/glew.h>
#include <GLFW/glfw3.h>

#include <png.h>

#include <fmt/core.h>

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

void MeshDiff::writeToPNG(const std::string & outPrefix)
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
        triangles.size() * sizeof(MeshDiff::GLTriangle),
        triangles.data(), GL_STATIC_DRAW);


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

    GLuint renderBuffer;
    glGenRenderbuffers(1, &renderBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, renderBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, width, width);
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, renderBuffer);
    if (glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
        throw std::runtime_error("Framebuffer is not complete.");
    glDisable(GL_DEPTH_TEST);

    glViewport(0, 0, width, width);


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
        std::cerr << fmt::format("Writing '{}'...\n", outFile);

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

    for (auto & tile : tiles) {

        glUniform1f(tileAttrib, tile);

        glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        std::cerr << fmt::format("Rendering tile {}...\n", tile);

        glDrawArrays(GL_TRIANGLES, 0, triangles.size() * 3);

        assert(glGetError() == GL_NO_ERROR);

        std::vector<Pixel> data(width * width);
        glReadBuffer(GL_COLOR_ATTACHMENT0);
        //glReadnPixels(0, 0, width, width, GL_RGB, GL_UNSIGNED_SHORT, data.size() * sizeof(Pixel), data.data());
        glReadPixels(0, 0, width, width, GL_RGB, GL_UNSIGNED_SHORT, data.data());

        if (auto err = glGetError())
            throw std::runtime_error("GL error: " + std::to_string(err));

        asyncs.push_back(std::async(std::launch::async, writePNG, tile, std::move(data)));
    }

    std::cerr << fmt::format("Waiting...\n");
    for (auto & async : asyncs)
        async.get();
}
