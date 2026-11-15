
#include "nonlinearfem/LoadStepping.h"
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
    std::cout << " Q4 Prescribed Displacement Test\n";
    std::cout << "====================================\n";

    // 1. Mesh

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

    // 2. Material

    MaterialParameters material{};

    material.youngModulus = 2.1e11;
    material.poissonRatio = 0.30;

    // 3. No external force.
    NodalLoads loads{};

    // 4. Boundary conditions.

    const double prescribedDisplacement = 1.0e-3;

    DisplacementConstraints constraints = {

        // Left boundary: ux = 0
        {0, 0, 0.0},
        {3, 0, 0.0},
        {6, 0, 0.0},

        // Remove rigid body motion.
        {0, 1, 0.0},

        // Right boundary: prescribed ux
        {2, 0, prescribedDisplacement},
        {5, 0, prescribedDisplacement},
        {8, 0, prescribedDisplacement}
    };

    GlobalDisplacement initialDisplacement(
        2 * mesh.nodes.size(), 0.0
    );

    // 5. Load stepping options.

    LoadSteppingOptions options{};

    options.initialIncrement = 0.25;
    options.minimumIncrement = 0.01;
    options.maximumIncrement = 0.25;

    options.cutbackFactor = 0.5;
    options.growthFactor = 1.0;

    options.maxRetries = 10;
    options.maxLoadSteps = 10;

    options.newtonOptions.maxIterations = 30;
    options.newtonOptions.residualTolerance = 1.0e-7;
    options.newtonOptions.displacementTolerance = 1.0e-10;

    // 6. Solve.

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

    // 7. Check convergence.

    std::cout << "\n========== Load History ==========\n";

    for (std::size_t i = 0; i < result.history.size(); ++i) {

        const auto& record = result.history[i];

        std::cout
            << "Step: " << i + 1
            << " | Lambda: " << record.loadFactor
            << " | Increment: " << record.loadIncrement
            << " | Retries: " << record.retries
            << " | Converged: " << record.converged
            << '\n';
    }

    // 8. Verify prescribed displacements.

    double constraintError = 0.0;

    for (const auto& constraint : constraints) {

        const std::size_t position =
            2 * constraint.node + constraint.direction;

        constraintError = std::max(
            constraintError,
            std::abs(
                result.displacement[position] -
                constraint.value
            ));
    }

    // 9. Verify final equilibrium.

    double relativeResidual = 0.0;
    double absoluteResidual = 0.0;

    if (result.converged) {

        const GlobalAssemblyResult assembly =
            assembleGlobalSystem(
                mesh,
                result.displacement,
                material,
                1.0
            );

        std::vector<bool> constrained(
            initialDisplacement.size(), false
        );

        for (const auto& constraint : constraints) {
            const std::size_t position =
                2 * constraint.node + constraint.direction;

            constrained[position] = true;
        }

        double residualSquared = 0.0;
        double internalForceSquared = 0.0;

        for (std::size_t i = 0;
            i < initialDisplacement.size();
             ++i) {

            if (!constrained[i]) {
                const double residual =
                    -assembly.internalForce[i];

                residualSquared += residual * residual;
            }

            // Include constrained DOFs as force reference.
            internalForceSquared +=
                assembly.internalForce[i] *
                assembly.internalForce[i];
        }

        absoluteResidual = std::sqrt(residualSquared);

        const double referenceForce =
            std::max(1.0, std::sqrt(internalForceSquared));

        relativeResidual =
            absoluteResidual / referenceForce;
    }

    // 10. Output.

    std::cout << "\n========== Verification ==========\n";

    std::cout << "Final load factor: "
              << result.finalLoadFactor << '\n';

    std::cout << "Node 2 ux: "
              << result.displacement[4] << '\n';

    std::cout << "Node 5 ux: "
              << result.displacement[10] << '\n';

    std::cout << "Node 8 ux: "
              << result.displacement[16] << '\n';

    std::cout << "Constraint error: "
              << constraintError << '\n';

    std::cout << "Free DOF residual norm: "
              << relativeResidual << '\n';

    std::cout << "Absolute residual: "
          << absoluteResidual << '\n';

    std::cout << "Relative residual: "
          << relativeResidual << '\n';

    const bool passed =
        result.converged &&
        std::abs(result.finalLoadFactor - 1.0) < 1.0e-12 &&
        constraintError < 1.0e-10 &&
        relativeResidual < 1.0e-7;

    std::cout << "\n====================================\n";
    std::cout << (passed ? "TEST PASSED" : "TEST FAILED");
    std::cout << "\n====================================\n";

    if (!passed) {
        std::cout << "Failure reason: "
                  << result.failureReason << '\n';
    }

    return passed ? 0 : 1;
}
