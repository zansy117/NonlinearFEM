
#include "nonlinearfem/BoundaryConditions.h"

#include <cmath>
#include <iostream>
#include <algorithm>

using namespace NonlinearFEM;

int main() {

    constexpr double tolerance = 1.0e-12;

    // --------------------------------
    // Test 1: External force assembly
    // --------------------------------

    const NodalLoads loads = {
        {0, 0, 100.0},
        {0, 1, 200.0},
        {1, 0, 300.0},
        {1, 0, 50.0}
    };

    const auto force = assembleExternalForce(2, loads);

    const GlobalInternalForce expectedForce = {
        100.0, 200.0,
        350.0, 0.0
    };

    double forceError = 0.0;

    for (std::size_t i = 0; i < 4; ++i) {
        forceError = std::max(
            forceError,
            std::abs(force[i] - expectedForce[i])
        );
    }

    // --------------------------------
    // Test 2: Displacement constraints
    // --------------------------------

    GlobalTangentStiffness K = {
        {4.0, 1.0, 2.0, 0.0},
        {1.0, 3.0, 0.0, 0.0},
        {2.0, 0.0, 5.0, 1.0},
        {0.0, 0.0, 1.0, 2.0}
    };

    GlobalInternalForce residual = {
        10.0, 20.0, 30.0, 40.0
    };

    const GlobalDisplacement displacement = {
        0.2, 0.0,
        0.0, 0.0
    };

    const DisplacementConstraints constraints = {
        {0, 0, 0.5}
    };

    applyDisplacementConstraints(
        K,
        residual,
        displacement,
        constraints
    );

    // Expected displacement increment = 0.5 - 0.2 = 0.3

    const GlobalInternalForce expectedResidual = {
        0.3,
        19.7,
        29.4,
        40.0
    };

    const GlobalTangentStiffness expectedK = {
        {1.0, 0.0, 0.0, 0.0},
        {0.0, 3.0, 0.0, 0.0},
        {0.0, 0.0, 5.0, 1.0},
        {0.0, 0.0, 1.0, 2.0}
    };

    double residualError = 0.0;
    double stiffnessError = 0.0;

    for (std::size_t i = 0; i < 4; ++i) {

        residualError = std::max(
            residualError,
            std::abs(residual[i] - expectedResidual[i])
        );

        for (std::size_t j = 0; j < 4; ++j) {

            stiffnessError = std::max(
                stiffnessError,
                std::abs(K[i][j] - expectedK[i][j])
            );
        }
    }

    // --------------------------------
    // Output
    // --------------------------------

    std::cout << "Boundary Conditions Test\n";
    std::cout << "------------------------------\n";

    std::cout << "Force error     : " << forceError << '\n';
    std::cout << "Residual error  : " << residualError << '\n';
    std::cout << "Stiffness error : " << stiffnessError << '\n';

    if (forceError < tolerance &&
        residualError < tolerance &&
        stiffnessError < tolerance) {

        std::cout << "TEST PASSED\n";
        return 0;
    }

    std::cout << "TEST FAILED\n";
    return 1;
}
