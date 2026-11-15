
#include "nonlinearfem/Solver.h"
#include "nonlinearfem/Material.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using namespace NonlinearFEM;

int main() {

    // ==================================
    // 1. Create Q4 mesh
    // ==================================

    Q4Mesh mesh{};

    mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {1.0, 1.0},
        {0.0, 1.0}
    };

    mesh.elements = {
        {0, 1, 2, 3}
    };

    // ==================================
    // 2. Material parameters
    // ==================================

    MaterialParameters material{};

    material.youngModulus = 2.1e11;
    material.poissonRatio = 0.30;

    const double thickness = 1.0;

    // ==================================
    // 3. External loads
    // ==================================

    const double totalForce = 1.0e8;

    NodalLoads loads = {
        {1, 0, totalForce / 2.0},
        {2, 0, totalForce / 2.0}
    };

    GlobalInternalForce externalForce =
        assembleExternalForce(mesh.nodes.size(), loads);

    // ==================================
    // 4. Boundary conditions
    // ==================================

    DisplacementConstraints constraints = {
        {0, 0, 0.0},
        {0, 1, 0.0},
        {3, 0, 0.0}
    };

    // ==================================
    // 5. Initial displacement
    // ==================================

    GlobalDisplacement initialDisplacement(
        2 * mesh.nodes.size(), 0.0
    );

    // ==================================
    // 6. Solver options
    // ==================================

    NewtonRaphsonOptions options{};

    options.maxIterations = 30;
    options.residualTolerance = 1.0e-8;
    options.displacementTolerance = 1.0e-10;

    // ==================================
    // 7. Solve
    // ==================================

    NewtonRaphsonResult result =
        solveNewtonRaphson(
            mesh,
            material,
            externalForce,
            constraints,
            initialDisplacement,
            options,
            thickness
        );

    std::cout << "========== Solver Result ==========\n";

    std::cout << "Converged: "
              << std::boolalpha
              << result.converged << '\n';

    std::cout << "Iterations: "
              << result.iterations << '\n';

    std::cout << "Residual norm: "
              << result.residualNorm << '\n';

    if (!result.converged) {

        std::cout << "Failure reason: "
                  << result.failureReason << '\n';

        return 1;
    }

    // ==================================
    // 8. Print nodal displacements
    // ==================================

    std::cout << "\n========== Nodal Displacements ==========\n";

    for (std::size_t i = 0; i < mesh.nodes.size(); ++i) {

        std::cout << "Node " << i
                  << ": ux = " << result.displacement[2 * i]
                  << ", uy = " << result.displacement[2 * i + 1]
                  << '\n';
    }

    // ==================================
    // 9. Verify displacement constraints
    // ==================================

    double constraintError = 0.0;

    for (const auto& constraint : constraints) {

        const std::size_t position =
            2 * constraint.node + constraint.direction;

        constraintError = std::max(
            constraintError,
            std::abs(
                result.displacement[position] - constraint.value
            )
        );
    }

    std::cout << "\nConstraint error: "
              << constraintError << '\n';

    if (constraintError > 1.0e-10) {
        std::cerr << "Constraint verification failed.\n";
        return 1;
    }

    // ==================================
    // 10. Independently verify equilibrium
    // ==================================

    GlobalAssemblyResult finalSystem =
        assembleGlobalSystem(
            mesh,
            result.displacement,
            material,
            thickness
        );

    double residualNormSquared = 0.0;
    double externalForceNormSquared = 0.0;

    std::vector<bool> constrained(
        externalForce.size(), false
    );

    for (const auto& constraint : constraints) {

        const std::size_t position =
            2 * constraint.node + constraint.direction;

        constrained[position] = true;
    }

    for (std::size_t i = 0; i < externalForce.size(); ++i) {

        if (constrained[i]) {
            continue;
        }

        const double residual =
            externalForce[i] - finalSystem.internalForce[i];

        residualNormSquared += residual * residual;

        externalForceNormSquared +=
            externalForce[i] * externalForce[i];
    }

    const double relativeResidual =
        std::sqrt(residualNormSquared) /
        std::max(1.0, std::sqrt(externalForceNormSquared));

    std::cout << "\nIndependent relative residual: "
              << relativeResidual << '\n';

    if (relativeResidual > 1.0e-7) {
        std::cerr << "Equilibrium verification failed.\n";
        return 1;
    }

    // ==================================
    // 11. Verify displacement direction
    // ==================================

    if (result.displacement[2 * 1] <= 0.0 ||
        result.displacement[2 * 2] <= 0.0) {

        std::cerr << "Displacement direction verification failed.\n";
        return 1;
    }

    std::cout << "\nAll solver tests passed.\n";

    return 0;
}
