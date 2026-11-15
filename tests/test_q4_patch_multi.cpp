
#include "nonlinearfem/LoadStepping.h"
#include "nonlinearfem/Assembly.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace NonlinearFEM;

int main() {

    std::cout << std::scientific
              << std::setprecision(10);

    std::cout << "====================================\n";
    std::cout << " Q4 Multi-Element Patch Test\n";
    std::cout << "====================================\n";

    // ==================================
    // 1. Mesh: 2 x 2 Q4 elements
    // ==================================

    Q4Mesh mesh{};

    mesh.nodes = {
        {0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0},
        {0.0, 0.5}, {0.5, 0.5}, {1.0, 0.5},
        {0.0, 1.0}, {0.5, 1.0}, {1.0, 1.0}
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

    const double thickness = 1.0;

    // ==================================
    // 3. Prescribed affine deformation
    // ==================================

    // F = [ 1.020   0.010 ]
    //     [ 0.005   0.990 ]

    const double F11 = 1.020;
    const double F12 = 0.010;
    const double F21 = 0.005;
    const double F22 = 0.990;

    // u = (F - I) X

    auto exactDisplacement = [&](
        const Q4NodeCoordinate& X
    ) -> std::array<double, 2> {

        const double x = X[0];
        const double y = X[1];

        return {
            (F11 - 1.0) * x + F12 * y,
            F21 * x + (F22 - 1.0) * y
        };
    };

    // ==================================
    // 4. Boundary constraints
    // ==================================

    // Node 4 is the center node.
    // All other nodes are boundary nodes.

    const std::array<std::size_t, 8> boundaryNodes = {
        0, 1, 2, 3, 5, 6, 7, 8
    };

    DisplacementConstraints constraints{};

    for (const std::size_t node : boundaryNodes) {

        const auto displacement =
            exactDisplacement(mesh.nodes[node]);

        constraints.push_back({
            node, 0, displacement[0]
        });

        constraints.push_back({
            node, 1, displacement[1]
        });
    }

    // ==================================
    // 5. External loads and initial state
    // ==================================

    NodalLoads loads{};

    GlobalDisplacement initialDisplacement(
        2 * mesh.nodes.size(), 0.0
    );

    // ==================================
    // 6. Load stepping options
    // ==================================

    LoadSteppingOptions options{};

    options.initialIncrement = 0.25;
    options.minimumIncrement = 0.01;
    options.maximumIncrement = 0.25;

    options.cutbackFactor = 0.5;
    options.growthFactor = 1.0;

    options.maxRetries = 10;
    options.maxLoadSteps = 10;

    options.newtonOptions.maxIterations = 30;
    options.newtonOptions.residualTolerance = 1.0e-9;
    options.newtonOptions.displacementTolerance = 1.0e-11;

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
            thickness
        );

    // ==================================
    // 8. Load history
    // ==================================

    std::cout << "\n========== Load History ==========\n";

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
    }

    // ==================================
    // 9. Check boundary displacement
    // ==================================

    double boundaryError = 0.0;

    for (const auto& constraint : constraints) {

        const std::size_t position =
            2 * constraint.node + constraint.direction;

        boundaryError = std::max(
            boundaryError,
            std::abs(
                result.displacement[position] -
                constraint.value
            )
        );
    }

    // ==================================
    // 10. Check center displacement
    // ==================================

    constexpr std::size_t centerNode = 4;

    const auto expectedCenter =
        exactDisplacement(mesh.nodes[centerNode]);

    const double numericalUx =
        result.displacement[2 * centerNode];

    const double numericalUy =
        result.displacement[2 * centerNode + 1];

    const double centerError = std::sqrt(
        std::pow(numericalUx - expectedCenter[0], 2) +
        std::pow(numericalUy - expectedCenter[1], 2)
    );

    // ==================================
    // 11. Independent equilibrium check
    // ==================================

    double absoluteResidual = 0.0;
    double relativeResidual = 0.0;

    if (result.converged) {

        const GlobalAssemblyResult assembly =
            assembleGlobalSystem(
                mesh,
                result.displacement,
                material,
                thickness
            );

        // Only the center node is free.

        double residualSquared = 0.0;
        double internalForceSquared = 0.0;

        for (std::size_t i = 0; i < assembly.internalForce.size(); ++i) {

            internalForceSquared +=
                assembly.internalForce[i] *
                assembly.internalForce[i];
        }

        for (std::size_t dof = 0; dof < 2; ++dof) {

            const std::size_t position =
                2 * centerNode + dof;

            residualSquared +=
                assembly.internalForce[position] *
                assembly.internalForce[position];
        }

        absoluteResidual = std::sqrt(residualSquared);

        relativeResidual =
            absoluteResidual /
            std::max(1.0, std::sqrt(internalForceSquared));
    }

    // ==================================
    // 12. Verification
    // ==================================

    std::cout << "\n========== Verification ==========\n";

    std::cout << "Final load factor: "
              << result.finalLoadFactor << '\n';

    std::cout << "\nCenter node displacement:\n";

    std::cout << "Numerical ux: "
              << numericalUx << '\n';

    std::cout << "Analytical ux: "
              << expectedCenter[0] << '\n';

    std::cout << "Numerical uy: "
              << numericalUy << '\n';

    std::cout << "Analytical uy: "
              << expectedCenter[1] << '\n';

    std::cout << "\nBoundary displacement error: "
              << boundaryError << '\n';

    std::cout << "Center displacement error: "
              << centerError << '\n';

    std::cout << "Center equilibrium residual: "
              << absoluteResidual << '\n';

    std::cout << "Relative residual: "
              << relativeResidual << '\n';

    const bool passed =
        result.converged &&
        std::abs(result.finalLoadFactor - 1.0) < 1.0e-12 &&
        boundaryError < 1.0e-10 &&
        centerError < 1.0e-9 &&
        relativeResidual < 1.0e-8;

    std::cout << "\n====================================\n";
    std::cout << (passed ? "TEST PASSED" : "TEST FAILED");
    std::cout << "\n====================================\n";

    if (!passed) {
        std::cout << "Failure reason: "
                  << result.failureReason << '\n';
    }

    return passed ? 0 : 1;
}
