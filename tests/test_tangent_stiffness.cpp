
#include "nonlinearfem/Assembly.h"
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
// Frobenius norm
// ==================================

double frobeniusNorm(
    const GlobalTangentStiffness& A
) {

    double sum = 0.0;

    for (const auto& row : A) {
        for (double value : row) {
            sum += value * value;
        }
    }

    return std::sqrt(sum);
}


// ==================================
// Relative tangent error
// ==================================

double relativeError(
    const GlobalTangentStiffness& analytical,
    const GlobalTangentStiffness& numerical
) {

    const std::size_t N = analytical.size();

    double errorSquared = 0.0;

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {

            const double difference =
                analytical[i][j] - numerical[i][j];

            errorSquared += difference * difference;
        }
    }

    return std::sqrt(errorSquared) /
        std::max(1.0, frobeniusNorm(numerical));
}


// ==================================
// Symmetry error
// ==================================

double symmetryError(
    const GlobalTangentStiffness& K
) {

    const std::size_t N = K.size();

    double errorSquared = 0.0;

    for (std::size_t i = 0; i < N; ++i) {
        for (std::size_t j = 0; j < N; ++j) {

            const double difference =
                K[i][j] - K[j][i];

            errorSquared += difference * difference;
        }
    }

    return std::sqrt(errorSquared) /
        std::max(1.0, frobeniusNorm(K));
}


// ==================================
// Central finite difference
// ==================================

GlobalTangentStiffness computeFiniteDifferenceTangent(
    const Q4Mesh& mesh,
    const MaterialParameters& material,
    const GlobalDisplacement& displacement,
    double thickness,
    double step
) {

    const std::size_t N = displacement.size();

    GlobalTangentStiffness Kfd(
        N, std::vector<double>(N, 0.0)
    );

    for (std::size_t j = 0; j < N; ++j) {

        GlobalDisplacement uPlus = displacement;
        GlobalDisplacement uMinus = displacement;

        uPlus[j] += step;
        uMinus[j] -= step;

        const GlobalAssemblyResult plus =
            assembleGlobalSystem(
                mesh, uPlus, material, thickness
            );

        const GlobalAssemblyResult minus =
            assembleGlobalSystem(
                mesh, uMinus, material, thickness
            );

        for (std::size_t i = 0; i < N; ++i) {

            Kfd[i][j] =
                (plus.internalForce[i] -
                 minus.internalForce[i])
                / (2.0 * step);
        }
    }

    return Kfd;
}


// ==================================
// Main test
// ==================================

int main() {

    std::cout << std::scientific
              << std::setprecision(8);

    std::cout
        << "========================================\n"
        << " Q4 Tangent Stiffness Verification\n"
        << "========================================\n";


    // 1. Mesh

    Q4Mesh mesh = create2x2Mesh();


    // 2. Material

    MaterialParameters material{};

    material.youngModulus = 2.1e11;
    material.poissonRatio = 0.30;

    const double thickness = 1.0;


    // 3. Nonzero displacement state

    GlobalDisplacement displacement = {

        0.0000,  0.0000,
        0.0100, -0.0020,
        0.0200, -0.0040,

        0.0020, -0.0050,
        0.0110, -0.0010,
        0.0220, -0.0030,

        0.0040, -0.0100,
        0.0130, -0.0060,
        0.0240, -0.0080
    };


    // 4. Analytical tangent stiffness

    const GlobalAssemblyResult analyticalSystem =
        assembleGlobalSystem(
            mesh,
            displacement,
            material,
            thickness
        );

    const GlobalTangentStiffness& analyticalK =
        analyticalSystem.tangentStiffness;


    // 5. Check analytical symmetry

    const double analyticalSymmetryError =
        symmetryError(analyticalK);

    std::cout << "\nAnalytical symmetry error : "
              << analyticalSymmetryError << '\n';


    // 6. Finite difference step sensitivity

    const std::vector<double> steps = {
        1.0e-5,
        3.0e-6,
        1.0e-6,
        3.0e-7,
        1.0e-7
    };

    std::cout
        << "\n========== Finite Difference Results ==========\n";

    std::cout
        << std::setw(16) << "Step"
        << std::setw(24) << "Relative Error"
        << std::setw(24) << "FD Symmetry Error"
        << '\n';

    double bestError = 1.0e100;
    double bestStep = 0.0;

    for (double step : steps) {

        const GlobalTangentStiffness numericalK =
            computeFiniteDifferenceTangent(
                mesh,
                material,
                displacement,
                thickness,
                step
            );

        const double error =
            relativeError(analyticalK, numericalK);

        const double fdSymmetry =
            symmetryError(numericalK);

        std::cout
            << std::setw(16) << step
            << std::setw(24) << error
            << std::setw(24) << fdSymmetry
            << '\n';

        if (error < bestError) {
            bestError = error;
            bestStep = step;
        }
    }


    // 7. Final result

    std::cout
        << "\n========================================\n";

    std::cout << "Best finite difference step : "
              << bestStep << '\n';

    std::cout << "Best relative error         : "
              << bestError << '\n';

    std::cout << "Analytical symmetry error   : "
              << analyticalSymmetryError << '\n';


    const bool passed =
        bestError < 1.0e-5 &&
        analyticalSymmetryError < 1.0e-10;

    if (passed) {
        std::cout << "\nTEST PASSED\n";
    }
    else {
        std::cout << "\nTEST FAILED\n";
    }

    std::cout
        << "========================================\n";

    return passed ? 0 : 1;
}
