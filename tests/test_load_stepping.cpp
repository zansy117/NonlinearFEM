#include "nonlinearfem/LoadStepping.h"
#include "nonlinearfem/BoundaryConditions.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <vector>

using namespace NonlinearFEM;

int main() {

    std::cout << std::scientific
              << std::setprecision(8);

    std::cout << "====================================\n";
    std::cout << " Q4 Load Stepping Test\n";
    std::cout << "====================================\n";


    // 1. Mesh

    Q4Mesh mesh{};

    mesh.nodes = {
        {0.0, 0.0},
        {0.5, 0.0},
        {1.0, 0.0},

        {0.0, 0.5},
        {0.5, 0.5},
        {1.0, 0.5},

        {0.0, 1.0},
        {0.5, 1.0},
        {1.0, 1.0}
    };

    mesh.elements = {
        {0, 1, 4, 3},
        {1, 2, 5, 4},
        {3, 4, 7, 6},
        {4, 5, 8, 7}
    };


    // 2. Material

    MaterialParameters material{};

    material.youngModulus = 2.1e11;
    material.poissonRatio = 0.30;


    // 3. Total external loads

    const double totalForce = 1.0e8;

    NodalLoads loads = {
        {2, 0, totalForce / 4.0},
        {5, 0, totalForce / 2.0},
        {8, 0, totalForce / 4.0}
    };


    // 4. Constraints

    DisplacementConstraints constraints = {
        {0, 0, 0.0},
        {3, 0, 0.0},
        {6, 0, 0.0},
        {0, 1, 0.0}
    };


    // 5. Initial displacement

    GlobalDisplacement initialDisplacement(
        2 * mesh.nodes.size(), 0.0
    );


    // 6. Load stepping options

    LoadSteppingOptions options{};

    options.initialIncrement = 0.25;
    options.minimumIncrement = 0.01;
    options.maximumIncrement = 0.25;

    options.cutbackFactor = 0.5;
    options.growthFactor = 1.5;

    options.maxRetries = 8;
    options.maxLoadSteps = 100;

    options.newtonOptions.maxIterations = 30;
    options.newtonOptions.residualTolerance = 1.0e-8;
    options.newtonOptions.displacementTolerance = 1.0e-10;


    // 7. Solve

    LoadSteppingResult result =
        solveLoadStepping(
            mesh,
            material,
            loads,
            constraints,
            initialDisplacement,
            options,
            1.0
        );


    // 8. Print load history

    std::cout << "\n========== Load History ==========\n";

    for (std::size_t i = 0;
         i < result.history.size();
         ++i) {

        const auto& record = result.history[i];

        std::cout << "Step: " << i + 1
                  << ", Lambda: " << record.loadFactor
                  << ", Increment: " << record.loadIncrement
                  << ", Iterations: " << record.iterations
                  << ", Retries: " << record.retries
                  << ", Residual: " << record.residualNorm
                  << ", Converged: " << record.converged
                  << '\n';
    }


    // 9. Final result

    std::cout << "\n========== Final Result ==========\n";

    std::cout << "Converged: "
              << result.converged << '\n';

    std::cout << "Final load factor: "
              << result.finalLoadFactor << '\n';

    if (!result.converged) {
        std::cout << "Failure: "
                  << result.failureReason << '\n';

        return 1;
    }


    // 10. Final displacement

    std::cout << "\n========== Nodal Displacements ==========\n";

    for (std::size_t i = 0;
         i < mesh.nodes.size();
         ++i) {

        std::cout << "Node " << i
                  << ": ux = "
                  << result.displacement[2 * i]
                  << ", uy = "
                  << result.displacement[2 * i + 1]
                  << '\n';
    }


    // 11. Verify final equilibrium

    const GlobalInternalForce finalExternalForce =
        assembleExternalForce(
            mesh.nodes.size(),
            loads
        );

    const GlobalAssemblyResult finalSystem =
        assembleGlobalSystem(
            mesh,
            result.displacement,
            material,
            1.0
        );

    std::vector<bool> constrained(
        finalExternalForce.size(), false
    );

    for (const auto& constraint : constraints) {

        const std::size_t position =
            2 * constraint.node + constraint.direction;

        constrained[position] = true;
    }

    double residualSquared = 0.0;
    double forceSquared = 0.0;

    for (std::size_t i = 0;
         i < finalExternalForce.size();
         ++i) {

        if (constrained[i]) continue;

        const double residual =
            finalExternalForce[i]
            - finalSystem.internalForce[i];

        residualSquared += residual * residual;
        forceSquared +=
            finalExternalForce[i] * finalExternalForce[i];
    }

    const double relativeResidual =
        std::sqrt(residualSquared) /
        std::max(1.0, std::sqrt(forceSquared));

    std::cout << "\nIndependent relative residual: "
              << relativeResidual << '\n';


    // 12. Verify final load factor

    const bool loadFactorPassed =
        std::abs(result.finalLoadFactor - 1.0) < 1.0e-12;

    const bool residualPassed =
        relativeResidual < 1.0e-7;

    const bool passed =
        result.converged &&
        loadFactorPassed &&
        residualPassed;

    std::cout << "\n====================================\n";

    std::cout << (passed ? "TEST PASSED" : "TEST FAILED");

    std::cout << "\n====================================\n";

    return passed ? 0 : 1;
}
