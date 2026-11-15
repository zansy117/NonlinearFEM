
#include "nonlinearfem/LoadStepping.h"
#include "nonlinearfem/BoundaryConditions.h"
#include "nonlinearfem/Assembly.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace NonlinearFEM;

int main() {

    std::cout << std::scientific
              << std::setprecision(8);

    std::cout << "====================================\n";
    std::cout << " Q4 Adaptive Cutback Test\n";
    std::cout << "====================================\n";


    // ==================================
    // 1. Mesh
    // ==================================

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


    // ==================================
    // 2. Material
    // ==================================

    MaterialParameters material{};

    material.youngModulus = 2.1e11;
    material.poissonRatio = 0.30;


    // ==================================
    // 3. External loads
    // ==================================

    const double totalForce = 1.0e11;

    NodalLoads loads = {
        {2, 0, totalForce / 4.0},
        {5, 0, totalForce / 2.0},
        {8, 0, totalForce / 4.0}
    };


    // ==================================
    // 4. Boundary conditions
    // ==================================

    DisplacementConstraints constraints = {
        {0, 0, 0.0},
        {3, 0, 0.0},
        {6, 0, 0.0},
        {0, 1, 0.0}
    };


    // ==================================
    // 5. Initial displacement
    // ==================================

    GlobalDisplacement initialDisplacement(
        2 * mesh.nodes.size(), 0.0
    );


    // ==================================
    // 6. Load stepping options
    // ==================================

    LoadSteppingOptions options{};

    options.initialIncrement = 1.0;
    options.minimumIncrement = 1.0e-4;
    options.maximumIncrement = 1.0;

    options.cutbackFactor = 0.5;
    options.growthFactor = 1.0;

    options.maxRetries = 15;
    options.maxLoadSteps = 2500;

    // Deliberately restrict Newton iterations.
    options.newtonOptions.maxIterations = 2;

    options.newtonOptions.residualTolerance = 1.0e-10;
    options.newtonOptions.displacementTolerance = 1.0e-10;


    // ==================================
    // 7. Solve
    // ==================================

    const LoadSteppingResult result =
        solveLoadStepping(
            mesh,
            material,
            loads,
            constraints,
            initialDisplacement,
            options,
            1.0
        );


    // ==================================
    // 8. Print history
    // ==================================

    std::cout << "\n========== Load History ==========\n";

    bool cutbackTriggered = false;

    for (std::size_t i = 0;
         i < result.history.size();
         ++i) {

        const auto& record = result.history[i];

        std::cout
            << "Step: " << i + 1
            << " | Lambda: " << record.loadFactor
            << " | Increment: " << record.loadIncrement
            << " | Retries: " << record.retries
            << " | Iterations: " << record.iterations
            << " | Residual: " << record.residualNorm
            << " | Converged: " << record.converged
            << '\n';

        if (record.retries > 0) {
            cutbackTriggered = true;
        }
    }


    // ==================================
    // 9. Verify final equilibrium
    // ==================================

    double relativeResidual = 0.0;

    if (result.converged) {

        const GlobalInternalForce externalForce =
            assembleExternalForce(
                mesh.nodes.size(),
                loads
            );

        const GlobalAssemblyResult assembly =
            assembleGlobalSystem(
                mesh,
                result.displacement,
                material,
                1.0
            );

        std::vector<bool> constrained(
            externalForce.size(), false
        );

        for (const auto& constraint : constraints) {

            const std::size_t position =
                2 * constraint.node + constraint.direction;

            constrained[position] = true;
        }

        double residualSquared = 0.0;
        double forceSquared = 0.0;

        for (std::size_t i = 0;
             i < externalForce.size();
             ++i) {

            if (constrained[i]) continue;

            const double residual =
                externalForce[i] -
                assembly.internalForce[i];

            residualSquared += residual * residual;

            forceSquared +=
                externalForce[i] * externalForce[i];
        }

        relativeResidual =
            std::sqrt(residualSquared) /
            std::max(1.0, std::sqrt(forceSquared));
    }


    // ==================================
    // 10. Final checks
    // ==================================

    const bool convergencePassed =
        result.converged;

    const bool loadFactorPassed =
        std::abs(result.finalLoadFactor - 1.0) < 1.0e-12;

    const bool cutbackPassed =
        cutbackTriggered;

    const bool residualPassed =
        relativeResidual < 1.0e-7;

    std::cout << "\n========== Verification ==========\n";

    std::cout << "Cutback triggered: "
              << cutbackTriggered << '\n';

    std::cout << "Final load factor: "
              << result.finalLoadFactor << '\n';

    std::cout << "Independent relative residual: "
              << relativeResidual << '\n';

    std::cout << "Converged: "
              << convergencePassed << '\n';

    const bool passed =
        convergencePassed &&
        loadFactorPassed &&
        cutbackPassed &&
        residualPassed;

    std::cout << "\n====================================\n";
    std::cout << (passed ? "TEST PASSED" : "TEST FAILED");
    std::cout << "\n====================================\n";

    if (!passed) {
        std::cout << "Failure reason: "
                  << result.failureReason << '\n';
    }

    return passed ? 0 : 1;
}
