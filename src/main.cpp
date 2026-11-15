#include "nonlinearfem/ModelIO.h"
#include <iostream>
#include <stdexcept>

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--help") {
        std::cout << "Usage: NonlinearFEM model.nfem output_directory\n"
                     "Output: results.html, nodes.csv, history.csv\n";
        return 0;
    }
    if (argc != 3) {
        std::cerr << "Usage: NonlinearFEM model.nfem output_directory\n";
        return 1;
    }
    try {
        const auto model = NonlinearFEM::readModelFile(argv[1]);
        if (std::filesystem::exists(argv[2]))
            throw std::runtime_error("Output directory already exists; choose a new run directory.");
        const auto result = NonlinearFEM::solveLoadStepping(
            model.mesh, model.material, model.loads, model.constraints,
            NonlinearFEM::GlobalDisplacement(2 * model.mesh.nodes.size(), 0.0),
            model.options, model.thickness);
        if (!result.converged) {
            std::cerr << "Analysis failed at load factor " << result.finalLoadFactor
                      << ": " << result.failureReason << '\n';
            return 2;
        }
        NonlinearFEM::writeResults(model, result, argv[2]);
        std::cout << "Converged: nodes=" << model.mesh.nodes.size()
                  << ", elements=" << model.mesh.elements.size()
                  << ", accepted steps=" << result.history.size()
                  << "\nOpen " << (std::filesystem::path(argv[2]) / "results.html").string()
                  << " in your browser.\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
