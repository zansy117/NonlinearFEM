
#include "nonlinearfem/ElementQ4.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>

using namespace NonlinearFEM;

namespace {

    constexpr double thickness = 1.0;

    constexpr double gauss = 0.57735026918962576451;

    using Vec2 = std::array<double, 2>;

    using ElementNodes = std::array<std::size_t, 4>;

    using Matrix8 = Q4TangentStiffness;

    const std::array<Vec2, 9> nodes = {{

        {0.00,  0.00},
        {0.55, -0.05},
        {1.15,  0.05},

        {-0.08, 0.48},
        {0.58,  0.57},
        {1.18,  0.52},

        {-0.02, 1.05},
        {0.62,  1.12},
        {1.25,  1.08}

    }};

    const std::array<ElementNodes, 4> elements = {{

        {0, 1, 4, 3},
        {1, 2, 5, 4},
        {3, 4, 7, 6},
        {4, 5, 8, 7}

    }};

    const std::array<double, 5> steps = {
        1.0e-4,
        1.0e-5,
        1.0e-6,
        1.0e-7,
        1.0e-8
    };

    double matrixNorm(const Matrix8& A) {

        double sum = 0.0;

        for (std::size_t i = 0; i < 8; ++i) {
            for (std::size_t j = 0; j < 8; ++j) {
                sum += A[i][j] * A[i][j];
            }
        }

        return std::sqrt(sum);
    }

    double relativeMatrixError(
        const Matrix8& A,
        const Matrix8& B
    ) {

        double errorSquared = 0.0;

        for (std::size_t i = 0; i < 8; ++i) {
            for (std::size_t j = 0; j < 8; ++j) {

                const double difference =
                    A[i][j] - B[i][j];

                errorSquared += difference * difference;
            }
        }

        return std::sqrt(errorSquared) /
               std::max(1.0, matrixNorm(A));
    }

    double symmetryError(const Matrix8& K) {

        double errorSquared = 0.0;

        for (std::size_t i = 0; i < 8; ++i) {
            for (std::size_t j = 0; j < 8; ++j) {

                const double difference =
                    K[i][j] - K[j][i];

                errorSquared += difference * difference;
            }
        }

        return std::sqrt(errorSquared) /
               std::max(1.0, matrixNorm(K));
    }

    Q4NodalValues getCoordinates(
        const ElementNodes& element
    ) {

        Q4NodalValues X{};

        for (std::size_t i = 0; i < 4; ++i) {
            X[i] = nodes[element[i]];
        }

        return X;
    }

    Q4NodalValues getDisplacements(
        const ElementNodes& element,
        const std::array<Vec2, 9>& globalDisplacement
    ) {

        Q4NodalValues u{};

        for (std::size_t i = 0; i < 4; ++i) {
            u[i] = globalDisplacement[element[i]];
        }

        return u;
    }

    Matrix8 computeFiniteDifferenceTangent(
        const Q4NodalValues& X,
        const Q4NodalValues& u,
        const MaterialParameters& material,
        double h
    ) {

        Matrix8 KFD{};

        for (std::size_t j = 0; j < 8; ++j) {

            Q4NodalValues uPlus = u;
            Q4NodalValues uMinus = u;

            const std::size_t node = j / 2;
            const std::size_t direction = j % 2;

            uPlus[node][direction] += h;
            uMinus[node][direction] -= h;

            const Q4InternalForce forcePlus =
                computeQ4InternalForce(
                    X, uPlus, material, thickness
                );

            const Q4InternalForce forceMinus =
                computeQ4InternalForce(
                    X, uMinus, material, thickness
                );

            for (std::size_t i = 0; i < 8; ++i) {

                KFD[i][j] =
                    (forcePlus[i] - forceMinus[i]) /
                    (2.0 * h);
            }
        }

        return KFD;
    }

}


int main() {

    std::cout << std::scientific
              << std::setprecision(10);

    std::cout
        << "====================================\n"
        << " Distorted Q4 Tangent Stiffness Test\n"
        << "====================================\n";


    // ==================================
    // 1. Material
    // ==================================

    MaterialParameters material{};

    material.youngModulus = 2.1e11;
    material.poissonRatio = 0.30;


    // ==================================
    // 2. Non-affine displacement field
    // ==================================

    const Matrix2D F = {{
        {1.020, 0.010},
        {0.005, 0.990}
    }};

    std::array<Vec2, 9> globalDisplacement{};

    for (std::size_t i = 0; i < nodes.size(); ++i) {

        const double x = nodes[i][0];
        const double y = nodes[i][1];

        // Affine deformation.

        globalDisplacement[i][0] =
            (F[0][0] - 1.0) * x +
            F[0][1] * y;

        globalDisplacement[i][1] =
            F[1][0] * x +
            (F[1][1] - 1.0) * y;

        // Add a small non-affine perturbation.

        globalDisplacement[i][0] +=
            1.0e-3 * std::sin(1.7 * x + 0.8 * y);

        globalDisplacement[i][1] +=
            8.0e-4 * std::cos(0.9 * x - 1.3 * y);
    }


    // ==================================
    // 3. Analytical symmetry verification
    // ==================================

    double maxSymmetryError = 0.0;

    for (std::size_t e = 0; e < elements.size(); ++e) {

        const Q4NodalValues X =
            getCoordinates(elements[e]);

        const Q4NodalValues u =
            getDisplacements(elements[e], globalDisplacement);

        const Matrix8 K =
            computeQ4TangentStiffness(
                X, u, material, thickness
            );

        maxSymmetryError = std::max(
            maxSymmetryError,
            symmetryError(K)
        );
    }


    // ==================================
    // 4. Finite difference verification
    // ==================================

    std::array<double, 5> errors{};

    for (std::size_t s = 0; s < steps.size(); ++s) {

        double maxElementError = 0.0;

        for (std::size_t e = 0; e < elements.size(); ++e) {

            const Q4NodalValues X =
                getCoordinates(elements[e]);

            const Q4NodalValues u =
                getDisplacements(elements[e], globalDisplacement);

            const Matrix8 analyticalK =
                computeQ4TangentStiffness(
                    X, u, material, thickness
                );

            const Matrix8 numericalK =
                computeFiniteDifferenceTangent(
                    X, u, material, steps[s]
                );

            const double error =
                relativeMatrixError(
                    analyticalK, numericalK
                );

            maxElementError = std::max(
                maxElementError, error
            );
        }

        errors[s] = maxElementError;

        std::cout
            << "Step size: " << steps[s]
            << " | Maximum relative error: "
            << errors[s] << '\n';
    }


    // ==================================
    // 5. Find best step
    // ==================================

    std::size_t bestIndex = 0;

    for (std::size_t i = 1; i < steps.size(); ++i) {

        if (errors[i] < errors[bestIndex]) {
            bestIndex = i;
        }
    }


    // ==================================
    // 6. Output
    // ==================================

    std::cout
        << "\n========== Verification ==========\n";

    std::cout
        << "Elements checked: "
        << elements.size() << '\n';

    std::cout
        << "DOFs per element: 8\n";

    std::cout
        << "Best finite difference step: "
        << steps[bestIndex] << '\n';

    std::cout
        << "Best relative error: "
        << errors[bestIndex] << '\n';

    std::cout
        << "Analytical symmetry error: "
        << maxSymmetryError << '\n';


    // ==================================
    // 7. Final checks
    // ==================================

    const bool passed =
        std::isfinite(errors[bestIndex]) &&
        errors[bestIndex] < 1.0e-6 &&
        maxSymmetryError < 1.0e-10;

    std::cout
        << "\n====================================\n"
        << (passed ? "TEST PASSED" : "TEST FAILED")
        << "\n====================================\n";

    return passed ? 0 : 1;
}
