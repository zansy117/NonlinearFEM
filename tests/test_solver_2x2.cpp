
#include "nonlinearfem/Solver.h"
#include "nonlinearfem/Material.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace NonlinearFEM;

// ==================================
// Create 2x2 Q4 mesh
// ==================================

Q4Mesh create2x2Mesh() {

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

    return mesh;
}


// ==================================
// Main test
// ==================================

int main() {

    std::cout << std::scientific
              << std::setprecision(8);

    std::cout << "====================================\n";
    std::cout << " 2x2 Q4 Nonlinear Solver Test\n";
    std::cout << "====================================\n\n";


    // ==================================
    // 1. Mesh
    // ==================================

    Q4Mesh mesh = create2x2Mesh();

    std::cout << "Node count    : "
              << mesh.nodes.size() << '\n';

    std::cout << "Element count : "
              << mesh.elements.size() << '\n';


    // ==================================
    // 2. Material
    // ==================================

    MaterialParameters material{};

    material.youngModulus = 2.1e11;
    material.poissonRatio = 0.30;

    const double thickness = 1.0;


    // ==================================
    // 3. External force
    // ==================================

    const double totalForce = 1.0e8;

    NodalLoads loads = {
        {2, 0, totalForce / 4.0},
        {5, 0, totalForce / 2.0},
        {8, 0, totalForce / 4.0}
    };

    GlobalInternalForce externalForce =
        assembleExternalForce(
            mesh.nodes.size(),
            loads
        );


    // ==================================
    // 4. Displacement constraints
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
    // 6. Newton-Raphson options
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


    std::cout << "\n========== Solver Result ==========\n";

    std::cout << "Converged      : "
              << std::boolalpha
              << result.converged << '\n';

    std::cout << "Iterations     : "
              << result.iterations << '\n';

    std::cout << "Residual norm  : "
              << result.residualNorm << '\n';


    if (!result.converged) {

        std::cout << "Failure reason : "
                  << result.failureReason << '\n';

        return 1;
    }


    // ==================================
    // 8. Nodal displacement
    // ==================================

    std::cout << "\n========== Nodal Displacements ==========\n";

    for (std::size_t i = 0;
         i < mesh.nodes.size();
         ++i) {

        std::cout << "Node " << i
                  << " : ux = "
                  << result.displacement[2 * i]
                  << ", uy = "
                  << result.displacement[2 * i + 1]
                  << '\n';
    }


    // ==================================
    // 9. Constraint verification
    // ==================================

    double constraintError = 0.0;

    std::vector<bool> constrained(
        externalForce.size(), false
    );

    for (const auto& constraint : constraints) {

        const std::size_t position =
            2 * constraint.node + constraint.direction;

        constrained[position] = true;

        constraintError = std::max(
            constraintError,
            std::abs(
                result.displacement[position]
                - constraint.value
            )
        );
    }

    std::cout << "\nConstraint error : "
              << constraintError << '\n';


    // ==================================
    // 10. Independent equilibrium test
    // ==================================

    GlobalAssemblyResult finalSystem =
        assembleGlobalSystem(
            mesh,
            result.displacement,
            material,
            thickness
        );

    double residualNormSquared = 0.0;
    double externalNormSquared = 0.0;

    for (std::size_t i = 0;
         i < externalForce.size();
         ++i) {

        if (constrained[i]) {
            continue;
        }

        const double residual =
            externalForce[i]
            - finalSystem.internalForce[i];

        residualNormSquared += residual * residual;

        externalNormSquared +=
            externalForce[i] * externalForce[i];
    }

    const double relativeResidual =
        std::sqrt(residualNormSquared)
        / std::max(
            1.0,
            std::sqrt(externalNormSquared)
        );

    std::cout << "Independent relative residual : "
              << relativeResidual << '\n';


    // ==================================
    // 11. Displacement direction
    // ==================================

    bool displacementDirectionPassed = true;

    for (std::size_t node : {2, 5, 8}) {

        if (result.displacement[2 * node] <= 0.0) {
            displacementDirectionPassed = false;
        }
    }

    std::cout << "Displacement direction : "
              << displacementDirectionPassed << '\n';


    // ==================================
    // 12. Global force balance
    // ==================================

    double supportReactionX = 0.0;
    double supportReactionY = 0.0;
    double externalForceX = 0.0;
    double externalForceY = 0.0;
    double appliedForceScale = 0.0;

    for (std::size_t node = 0; node < mesh.nodes.size(); ++node) {
        const std::size_t x = 2 * node;
        const std::size_t y = x + 1;

        externalForceX += externalForce[x];
        externalForceY += externalForce[y];
        appliedForceScale += std::abs(externalForce[x])
                           + std::abs(externalForce[y]);

        // Reactions exist only at constrained DOFs.
        if (constrained[x]) {
            supportReactionX +=
                finalSystem.internalForce[x] - externalForce[x];
        }
        if (constrained[y]) {
            supportReactionY +=
                finalSystem.internalForce[y] - externalForce[y];
        }
    }

    // Use the applied load scale for both directions, including zero Y load.
    const double forceScale = std::max(1.0, appliedForceScale);
    const double balanceErrorX =
        std::abs(supportReactionX + externalForceX) / forceScale;
    const double balanceErrorY =
        std::abs(supportReactionY + externalForceY) / forceScale;
    const bool forceBalancePassed =
        std::isfinite(balanceErrorX) && std::isfinite(balanceErrorY) &&
        balanceErrorX < 1.0e-7 && balanceErrorY < 1.0e-7;

    std::cout << "\n========== Global Force Balance ==========\n";
    std::cout << "External force X : " << externalForceX << '\n';
    std::cout << "Support reaction X : " << supportReactionX << '\n';
    std::cout << "External force Y : " << externalForceY << '\n';
    std::cout << "Support reaction Y : " << supportReactionY << '\n';
    std::cout << "Relative balance error X : " << balanceErrorX << '\n';
    std::cout << "Relative balance error Y : " << balanceErrorY << '\n';
    std::cout << "Force balance passed : " << forceBalancePassed << '\n';

    // ==================================
    // 13. Final verification
    // ==================================

    const bool passed =
        result.converged &&
        constraintError < 1.0e-10 &&
        relativeResidual < 1.0e-7 &&
        displacementDirectionPassed &&
        forceBalancePassed;

    std::cout << "\n====================================\n";

    if (passed) {
        std::cout << "TEST PASSED\n";
    }
    else {
        std::cout << "TEST FAILED\n";
    }

    std::cout << "====================================\n";

    return passed ? 0 : 1;
}
