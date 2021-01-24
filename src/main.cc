#include "mesh.hh"
#include "dhdm.hh"
#include "diff.hh"

#include <vector>
#include <string>
#include <iostream>
#include <future>

#include <fmt/core.h>

unsigned int applyMorphs(Mesh & mesh, size_t nrPaths, char * paths[])
{
    std::vector<std::pair<double, Dhdm>> dhdms;
    unsigned int level = 1;

    for (size_t i = 0; i < nrPaths; ++i) {
        std::string fn = paths[i];
        double weight = 1.0;
        auto c = fn.find('=');
        if (c != fn.npos) {
            weight = stod(fn.substr(c + 1));
            fn = std::string(fn, 0, c);
        }

        if (fn.size() >= 4 && std::string(fn, fn.size() - 4) == ".dsf") {
            if (auto dhdm = mesh.applyMorph(weight, fn))
                dhdms.push_back({weight, Dhdm(*dhdm)});
        } else {
            auto dhdm = Dhdm(fn);
            level = std::max(level, (unsigned int) dhdm.levels.size());
            dhdms.push_back({weight, std::move(dhdm)});
        }
    }

    mesh.subdivide(level, dhdms);

    return level;
}

int main(int argc, char * * argv)
{
    try {

        if (argc < 2) {
            std::cerr <<
                "Syntax: dhdm morphs-to-displacement <prefix> <base-mesh>.dsf <uv-map>.dsf [<morph>.[dsf|dhdm][=strength]]...\n"
                "           | morphs-to-obj <base-mesh>.dsf <uv-map>.dsf [<morph>.[dsf|dhdm][=strength]]... > dest.obj\n"
                "           | diff-objs <prefix> <base-mesh>.obj [<mesh>.obj[=strength]]...\n";
            return 1;
        }

        std::string verb = argv[1];

        if (verb == "morphs-to-displacement") {
            if (argc < 5)
                throw std::runtime_error("Missing arguments.");

            std::string prefix = argv[2];
            std::string baseMeshPath = argv[3];
            std::string uvMapPath = argv[4];

            auto baseMesh = Mesh::fromDSF(baseMeshPath, uvMapPath);
            auto finalMesh = baseMesh;
            auto level = applyMorphs(finalMesh, argc - 5, argv + 5);

            baseMesh.subdivide(level, {});

            std::cerr << fmt::format("Triangulating...\n");
            baseMesh.triangulate();
            finalMesh.triangulate();

            auto diff = MeshDiff::create(1.0, baseMesh, {{1.0, std::move(finalMesh)}});

            diff.writeToPNG(prefix);
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

            auto diff = MeshDiff::create(1.0, baseMesh, morphMeshes);

            diff.writeToPNG(prefix);
        }

        else
            throw std::runtime_error(fmt::format("Unknown command '{}'.", verb));

        return 0;
    } catch (std::exception & e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
