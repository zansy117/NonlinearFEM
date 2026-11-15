
#include "nonlinearfem/ElementQ4.h"
#include "nonlinearfem/Kinematics2D.h"
#include "nonlinearfem/Material.h"
#include "nonlinearfem/BoundaryConditions.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace NonlinearFEM;

namespace {

    constexpr double tolerance = 1.0e-9;
    constexpr double gauss = 0.57735026918962576451;
    constexpr double thickness = 1.0;

    using Vec2 = std::array<double, 2>;

    using Matrix = Matrix2D;

    using ElementNodes = std::array<std::size_t, 4>;

    const std::array<double, 2> gaussPoints = {
        -gauss, gauss
    };

    // ------------------------------------------------
    // Exact deformation gradient
    // ------------------------------------------------

    const Matrix exactF = {{
        {1.020, 0.010},
        {0.005, 0.990}
    }};

    // ------------------------------------------------
    // Error utilities
    // ------------------------------------------------

    double relativeError(double numerical, double analytical) {

        return std::abs(numerical - analytical) /
               std::max(1.0, std::abs(analytical));
    }

    // ------------------------------------------------
    // Analytical displacement
    // ------------------------------------------------

    Vec2 exactDisplacement(const Vec2& X) {

        return {
            (exactF[0][0] - 1.0) * X[0]
                + exactF[0][1] * X[1],

            exactF[1][0] * X[0]
                + (exactF[1][1] - 1.0) * X[1]
        };
    }

    // ------------------------------------------------
    // Independent shape function derivatives
    // ------------------------------------------------

    void computeShapeDerivatives(
        double xi,
        double eta,
        double hx,
        double hy,
        std::array<double, 4>& dNdx,
        std::array<double, 4>& dNdy
    ) {

        const std::array<double, 4> dNdxi = {
            (eta - 1.0) / 4.0,
            (1.0 - eta) / 4.0,
            (eta + 1.0) / 4.0,
            (-eta - 1.0) / 4.0
        };

        const std::array<double, 4> dNdeta = {
            (xi - 1.0) / 4.0,
            (-1.0 - xi) / 4.0,
            (xi + 1.0) / 4.0,
            (1.0 - xi) / 4.0
        };

        for (std::size_t i = 0; i < 4; ++i) {

            dNdx[i] = dNdxi[i] * 2.0 / hx;

            dNdy[i] = dNdeta[i] * 2.0 / hy;
        }
    }

    // ------------------------------------------------
    // Independent displacement gradient
    // ------------------------------------------------

    Matrix computeIndependentF(
        const Q4NodalValues& displacement,
        const std::array<double, 4>& dNdx,
        const std::array<double, 4>& dNdy
    ) {

        Matrix gradient{};

        for (std::size_t i = 0; i < 4; ++i) {

            gradient[0][0] += dNdx[i] * displacement[i][0];
            gradient[0][1] += dNdy[i] * displacement[i][0];

            gradient[1][0] += dNdx[i] * displacement[i][1];
            gradient[1][1] += dNdy[i] * displacement[i][1];
        }

        Matrix F = gradient;

        F[0][0] += 1.0;
        F[1][1] += 1.0;

        return F;
    }

    // ------------------------------------------------
    // Analytical Green-Lagrange strain
    // ------------------------------------------------

    const std::array<double, 3> exactStrain = {
        0.0202125,
        -0.0099,
        0.01515
    };

    // ------------------------------------------------
    // Analytical PK2 stress
    // ------------------------------------------------

    std::array<double, 3> computeExactStress(
        const MaterialParameters& material
    ) {

        const double E = material.youngModulus;
        const double nu = material.poissonRatio;

        return {
            E * (exactStrain[0] + nu * exactStrain[1])
                / (1.0 - nu * nu),

            E * (exactStrain[1] + nu * exactStrain[0])
                / (1.0 - nu * nu),

            E * exactStrain[2]
                / (2.0 * (1.0 + nu))
        };
    }

    // ------------------------------------------------
    // Analytical PK1 stress: P = F * S
    // ------------------------------------------------

    Matrix computeExactPK1(
        const std::array<double, 3>& S
    ) {

        Matrix P{};

        P[0][0] = exactF[0][0] * S[0]
                 + exactF[0][1] * S[2];

        P[0][1] = exactF[0][0] * S[2]
                 + exactF[0][1] * S[1];

        P[1][0] = exactF[1][0] * S[0]
                 + exactF[1][1] * S[2];

        P[1][1] = exactF[1][0] * S[2]
                 + exactF[1][1] * S[1];

        return P;
    }

    // ------------------------------------------------
    // Independent internal force integration
    // ------------------------------------------------

    Q4InternalForce computeIndependentInternalForce(
        const Q4NodalValues& coordinates,
        const Q4NodalValues& displacement,
        const Matrix& P
    ) {

        Q4InternalForce force{};

        const double hx =
            coordinates[1][0] - coordinates[0][0];

        const double hy =
            coordinates[3][1] - coordinates[0][1];

        const double detJ = hx * hy / 4.0;

        for (double xi : gaussPoints) {

            for (double eta : gaussPoints) {

                std::array<double, 4> dNdx{};
                std::array<double, 4> dNdy{};

                computeShapeDerivatives(
                    xi, eta, hx, hy,
                    dNdx, dNdy
                );

                const double weight = thickness * detJ;

                for (std::size_t i = 0; i < 4; ++i) {

                    force[2 * i] +=
                        (dNdx[i] * P[0][0]
                       + dNdy[i] * P[0][1]) * weight;

                    force[2 * i + 1] +=
                        (dNdx[i] * P[1][0]
                       + dNdy[i] * P[1][1]) * weight;
                }
            }
        }

        return force;
    }

}


// ====================================================
// Main
// ====================================================

int main() {

    std::cout << std::scientific
              << std::setprecision(10);

    std::cout
        << "====================================\n"
        << " Q4 Multi-Element Stress Patch Test\n"
        << "====================================\n";


    // ------------------------------------------------
    // 1. Mesh
    // ------------------------------------------------

    const std::array<Vec2, 9> nodes = {{

        {0.0, 0.0}, {0.5, 0.0}, {1.0, 0.0},

        {0.0, 0.5}, {0.5, 0.5}, {1.0, 0.5},

        {0.0, 1.0}, {0.5, 1.0}, {1.0, 1.0}

    }};

    const std::array<ElementNodes, 4> elements = {{

        {0, 1, 4, 3},
        {1, 2, 5, 4},
        {3, 4, 7, 6},
        {4, 5, 8, 7}

    }};


    // ------------------------------------------------
    // 2. Material
    // ------------------------------------------------

    MaterialParameters material{};

    material.youngModulus = 2.1e11;
    material.poissonRatio = 0.30;

    const auto exactS = computeExactStress(material);

    const Matrix exactP = computeExactPK1(exactS);


    // ------------------------------------------------
    // 3. Global force accumulation
    // ------------------------------------------------

    GlobalDisplacement dummy;

    std::array<double, 18> globalInternalForce{};


    // ------------------------------------------------
    // 4. Error statistics
    // ------------------------------------------------

    double maxFError = 0.0;
    double maxStrainError = 0.0;
    double maxStressError = 0.0;
    double maxInternalForceError = 0.0;

    std::size_t gaussPointCount = 0;


    // ------------------------------------------------
    // 5. Element loop
    // ------------------------------------------------

    for (std::size_t e = 0; e < elements.size(); ++e) {

        Q4NodalValues coordinates{};
        Q4NodalValues displacement{};

        for (std::size_t i = 0; i < 4; ++i) {

            const std::size_t node = elements[e][i];

            coordinates[i] = nodes[node];

            displacement[i] =
                exactDisplacement(nodes[node]);
        }

        const double hx =
            coordinates[1][0] - coordinates[0][0];

        const double hy =
            coordinates[3][1] - coordinates[0][1];


        // --------------------------------------------
        // 5.1 Gauss point verification
        // --------------------------------------------

        for (double xi : gaussPoints) {

            for (double eta : gaussPoints) {

                std::array<double, 4> dNdx{};
                std::array<double, 4> dNdy{};

                computeShapeDerivatives(
                    xi, eta, hx, hy,
                    dNdx, dNdy
                );

                const Matrix numericalF =
                    computeIndependentF(
                        displacement, dNdx, dNdy
                    );

                for (std::size_t i = 0; i < 2; ++i) {

                    for (std::size_t j = 0; j < 2; ++j) {

                        maxFError = std::max(
                            maxFError,
                            relativeError(
                                numericalF[i][j],
                                exactF[i][j]
                            )
                        );
                    }
                }

                const GreenLagrangeStrain strain =
                    computeGreenLagrangeStrain(numericalF);

                maxStrainError = std::max({

                    maxStrainError,

                    relativeError(strain.exx, exactStrain[0]),

                    relativeError(strain.eyy, exactStrain[1]),

                    relativeError(strain.exy, exactStrain[2])

                });

                const SecondPKStress stress =
                    computeSecondPKStress(strain, material);

                maxStressError = std::max({

                    maxStressError,

                    relativeError(stress.sxx, exactS[0]),

                    relativeError(stress.syy, exactS[1]),

                    relativeError(stress.sxy, exactS[2])

                });

                ++gaussPointCount;
            }
        }


        // --------------------------------------------
        // 5.2 Production internal force
        // --------------------------------------------

        const Q4InternalForce numericalForce =
            computeQ4InternalForce(
                coordinates,
                displacement,
                material,
                thickness
            );


        // --------------------------------------------
        // 5.3 Independent internal force
        // --------------------------------------------

        const Q4InternalForce analyticalForce =
            computeIndependentInternalForce(
                coordinates,
                displacement,
                exactP
            );

        double elementForceScale = 1.0;

        for (std::size_t i = 0; i < 8; ++i) {

            elementForceScale = std::max(
                elementForceScale,
                std::abs(analyticalForce[i])
            );
        }

        for (std::size_t i = 0; i < 8; ++i) {

            maxInternalForceError = std::max(
                maxInternalForceError,
                std::abs(
                    numericalForce[i] - analyticalForce[i]
                ) / elementForceScale
            );
        }


        // --------------------------------------------
        // 5.4 Assemble global internal force
        // --------------------------------------------

        for (std::size_t i = 0; i < 4; ++i) {

            const std::size_t node = elements[e][i];

            globalInternalForce[2 * node] +=
                numericalForce[2 * i];

            globalInternalForce[2 * node + 1] +=
                numericalForce[2 * i + 1];
        }

        std::cout
            << "Element " << e
            << " verified.\n";
    }


    // ------------------------------------------------
    // 6. Center node equilibrium
    // ------------------------------------------------

    const double centerResidual = std::hypot(
        globalInternalForce[2 * 4],
        globalInternalForce[2 * 4 + 1]
    );

    double globalForceScale = 1.0;

    for (double force : globalInternalForce) {

        globalForceScale = std::max(
            globalForceScale,
            std::abs(force)
        );
    }

    const double relativeCenterResidual =
        centerResidual / globalForceScale;


    // ------------------------------------------------
    // 7. Output
    // ------------------------------------------------

    std::cout
        << "\n========== Analytical Reference ==========\n";

    std::cout << "F:\n";

    for (const auto& row : exactF) {

        std::cout << row[0] << " "
                  << row[1] << '\n';
    }

    std::cout << "\nGreen-Lagrange strain:\n";

    for (double value : exactStrain) {

        std::cout << value << '\n';
    }

    std::cout << "\nPK2 stress:\n";

    for (double value : exactS) {

        std::cout << value << '\n';
    }


    std::cout
        << "\n========== Verification ==========\n";

    std::cout
        << "Gauss points checked: "
        << gaussPointCount << '\n';

    std::cout
        << "Maximum F relative error: "
        << maxFError << '\n';

    std::cout
        << "Maximum strain relative error: "
        << maxStrainError << '\n';

    std::cout
        << "Maximum PK2 stress relative error: "
        << maxStressError << '\n';

    std::cout
        << "Maximum element force relative error: "
        << maxInternalForceError << '\n';

    std::cout
        << "Center absolute residual: "
        << centerResidual << '\n';

    std::cout
        << "Center relative residual: "
        << relativeCenterResidual << '\n';


    // ------------------------------------------------
    // 8. Final checks
    // ------------------------------------------------

    const bool passed =

        maxFError < tolerance &&

        maxStrainError < tolerance &&

        maxStressError < tolerance &&

        maxInternalForceError < tolerance &&

        relativeCenterResidual < tolerance;


    std::cout
        << "\n====================================\n"
        << (passed ? "TEST PASSED" : "TEST FAILED")
        << "\n====================================\n";

    return passed ? 0 : 1;
}
