
#include "nonlinearfem/ElementQ4.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

using namespace NonlinearFEM;

int main() {

    // ---------------------------------
    // 1. Reference geometry
    // ---------------------------------

    const Q4NodalValues coordinates{{
        {{0.0, 0.0}},
        {{2.0, 0.0}},
        {{2.0, 1.0}},
        {{0.0, 1.0}}
    }};

    // ---------------------------------
    // 2. Material parameters
    // ---------------------------------

    const MaterialParameters material{
        2.1e5,
        0.30
    };

    const double thickness = 1.0;

    // ---------------------------------
    // 3. Non-zero displacement state
    // ---------------------------------

    const Q4NodalValues displacement{{
        {{0.001,  0.002}},
        {{0.015, -0.003}},
        {{0.020,  0.012}},
        {{-0.002, 0.008}}
    }};

    // ---------------------------------
    // 4. Analytical tangent stiffness
    // ---------------------------------

    const Q4TangentStiffness K =
        computeQ4TangentStiffness(
            coordinates,
            displacement,
            material,
            thickness
        );

    // ---------------------------------
    // 5. Finite difference tangent
    // ---------------------------------

    std::array<
        std::array<double, 8>,
        8
    > Kfd{};

    constexpr double h = 1.0e-7;

    for (std::size_t j = 0; j < 8; ++j) {

        Q4NodalValues uPlus = displacement;
        Q4NodalValues uMinus = displacement;

        const std::size_t node = j / 2;
        const std::size_t component = j % 2;

        uPlus[node][component] += h;
        uMinus[node][component] -= h;

        const Q4InternalForce fPlus =
            computeQ4InternalForce(
                coordinates,
                uPlus,
                material,
                thickness
            );

        const Q4InternalForce fMinus =
            computeQ4InternalForce(
                coordinates,
                uMinus,
                material,
                thickness
            );

        for (std::size_t i = 0; i < 8; ++i) {
            Kfd[i][j] =
                (fPlus[i] - fMinus[i]) / (2.0 * h);
        }
    }

    // ---------------------------------
    // 6. Relative Frobenius error
    // ---------------------------------

    double errorSquared = 0.0;
    double referenceSquared = 0.0;
    double symmetryErrorSquared = 0.0;

    for (std::size_t i = 0; i < 8; ++i) {
        for (std::size_t j = 0; j < 8; ++j) {

            const double difference = K[i][j] - Kfd[i][j];

            errorSquared += difference * difference;
            referenceSquared += Kfd[i][j] * Kfd[i][j];

            const double asymmetry = K[i][j] - K[j][i];

            symmetryErrorSquared += asymmetry * asymmetry;
        }
    }

    const double relativeError =
        std::sqrt(errorSquared) /
        std::max(1.0, std::sqrt(referenceSquared));

    const double symmetryError =
        std::sqrt(symmetryErrorSquared) /
        std::max(1.0, std::sqrt(referenceSquared));

    // ---------------------------------
    // 7. Output
    // ---------------------------------

    std::cout << std::scientific << std::setprecision(8);

    std::cout << "Finite Difference Tangent Test\n";
    std::cout << "--------------------------------\n";

    std::cout << "Relative error : "
              << relativeError << '\n';

    std::cout << "Symmetry error : "
              << symmetryError << '\n';

    constexpr double tolerance = 1.0e-6;

    if (relativeError < tolerance &&
        symmetryError < tolerance) {

        std::cout << "TEST PASSED\n";
        return 0;
    }

    std::cout << "TEST FAILED\n";
    return 1;
}
