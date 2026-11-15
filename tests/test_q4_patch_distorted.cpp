
#include "nonlinearfem/LoadStepping.h"
#include "nonlinearfem/Assembly.h"
#include "nonlinearfem/ElementQ4.h"
#include "nonlinearfem/Kinematics2D.h"
#include "nonlinearfem/Material.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace NonlinearFEM;

namespace {

    constexpr double g = 0.57735026918962576451;
    constexpr double thickness = 1.0;
    constexpr double tolerance = 1.0e-9;

    using Vec2 = std::array<double, 2>;
    using Matrix = Matrix2D;

    const std::array<double, 2> gaussPoints = {-g, g};

    const Matrix exactF = {{
        {1.020, 0.010},
        {0.005, 0.990}
    }};

    const std::array<double, 3> exactStrain = {
        0.0202125,
        -0.0099,
        0.01515
    };

    double relativeError(double a, double b) {
        return std::abs(a - b) /
               std::max(1.0, std::abs(b));
    }

    Vec2 exactDisplacement(const Vec2& X) {
        return {
            (exactF[0][0] - 1.0) * X[0]
                + exactF[0][1] * X[1],

            exactF[1][0] * X[0]
                + (exactF[1][1] - 1.0) * X[1]
        };
    }

    struct ShapeData {
        std::array<double, 4> dNdxi{};
        std::array<double, 4> dNdeta{};
        std::array<double, 4> dNdx{};
        std::array<double, 4> dNdy{};
        double detJ = 0.0;
    };

    ShapeData computeShapeData(
        const Q4NodalValues& X,
        double xi,
        double eta
    ) {

        ShapeData data{};

        data.dNdxi = {
            (eta - 1.0) / 4.0,
            (1.0 - eta) / 4.0,
            (1.0 + eta) / 4.0,
            (-1.0 - eta) / 4.0
        };

        data.dNdeta = {
            (xi - 1.0) / 4.0,
            (-1.0 - xi) / 4.0,
            (1.0 + xi) / 4.0,
            (1.0 - xi) / 4.0
        };

        double J00 = 0.0;
        double J01 = 0.0;
        double J10 = 0.0;
        double J11 = 0.0;

        for (std::size_t i = 0; i < 4; ++i) {

            J00 += data.dNdxi[i] * X[i][0];
            J01 += data.dNdeta[i] * X[i][0];

            J10 += data.dNdxi[i] * X[i][1];
            J11 += data.dNdeta[i] * X[i][1];
        }

        data.detJ = J00 * J11 - J01 * J10;

        if (data.detJ <= 0.0) {
            throw std::runtime_error(
                "Non-positive Jacobian determinant."
            );
        }

        // grad_X(N) = inverse(J)^T * grad_xi(N)

        for (std::size_t i = 0; i < 4; ++i) {

            data.dNdx[i] =
                (J11 * data.dNdxi[i]
                - J10 * data.dNdeta[i]) / data.detJ;

            data.dNdy[i] =
                (-J01 * data.dNdxi[i]
                + J00 * data.dNdeta[i]) / data.detJ;
        }

        return data;
    }

    Matrix computeIndependentF(
        const Q4NodalValues& displacement,
        const ShapeData& shape
    ) {

        Matrix F{};

        F[0][0] = 1.0;
        F[1][1] = 1.0;

        for (std::size_t i = 0; i < 4; ++i) {

            F[0][0] += shape.dNdx[i] * displacement[i][0];
            F[0][1] += shape.dNdy[i] * displacement[i][0];

            F[1][0] += shape.dNdx[i] * displacement[i][1];
            F[1][1] += shape.dNdy[i] * displacement[i][1];
        }

        return F;
    }

    std::array<double, 3> exactPK2(
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

    Matrix exactPK1(
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

    Q4InternalForce independentInternalForce(
        const Q4NodalValues& X,
        const Matrix& P
    ) {

        Q4InternalForce force{};

        for (double xi : gaussPoints) {

            for (double eta : gaussPoints) {

                const ShapeData shape =
                    computeShapeData(X, xi, eta);

                const double weight =
                    thickness * shape.detJ;

                for (std::size_t i = 0; i < 4; ++i) {

                    force[2 * i] +=
                        (shape.dNdx[i] * P[0][0]
                       + shape.dNdy[i] * P[0][1]) * weight;

                    force[2 * i + 1] +=
                        (shape.dNdx[i] * P[1][0]
                       + shape.dNdy[i] * P[1][1]) * weight;
                }
            }
        }

        return force;
    }

}


int main() {

    std::cout << std::scientific
              << std::setprecision(10);

    std::cout
        << "====================================\n"
        << " Distorted Q4 Multi-Element Patch Test\n"
        << "====================================\n";

    // ==================================
    // 1. Distorted mesh
    // ==================================

    Q4Mesh mesh{};

    mesh.nodes = {
        {0.00,  0.00},
        {0.55, -0.05},
        {1.15,  0.05},

        {-0.08, 0.48},
        {0.58,  0.57},
        {1.18,  0.52},

        {-0.02, 1.05},
        {0.62,  1.12},
        {1.25,  1.08}
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
    // 3. Prescribed affine displacement
    // ==================================

    const std::array<std::size_t, 8> boundaryNodes = {
        0, 1, 2, 3, 5, 6, 7, 8
    };

    DisplacementConstraints constraints{};

    for (const std::size_t node : boundaryNodes) {

        const Vec2 X = mesh.nodes[node];

        const Vec2 u = exactDisplacement(X);

        constraints.push_back({node, 0, u[0]});
        constraints.push_back({node, 1, u[1]});
    }

    NodalLoads loads{};

    GlobalDisplacement initialDisplacement(
        2 * mesh.nodes.size(), 0.0
    );

    // ==================================
    // 4. Load stepping
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

    if (!result.converged) {

        std::cout << "Solver failed: "
                  << result.failureReason << '\n';

        return 1;
    }

    // ==================================
    // 5. Boundary and center displacement
    // ==================================

    double boundaryError = 0.0;

    for (const auto& constraint : constraints) {

        const std::size_t dof =
            2 * constraint.node + constraint.direction;

        boundaryError = std::max(
            boundaryError,
            std::abs(
                result.displacement[dof] - constraint.value
            )
        );
    }

    const Vec2 expectedCenter =
        exactDisplacement(mesh.nodes[4]);

    const double centerUx = result.displacement[8];
    const double centerUy = result.displacement[9];

    const double centerError = std::hypot(
        centerUx - expectedCenter[0],
        centerUy - expectedCenter[1]
    );

    // ==================================
    // 6. Analytical stresses
    // ==================================

    const auto S = exactPK2(material);
    const Matrix P = exactPK1(S);

    // ==================================
    // 7. Element and Gauss point checks
    // ==================================

    double maxJacobianError = 0.0;
    double minDetJ = 1.0e100;

    double maxFError = 0.0;
    double maxStrainError = 0.0;
    double maxStressError = 0.0;
    double maxForceError = 0.0;

    std::array<double, 18> globalInternalForce{};

    std::size_t gaussCount = 0;

    for (std::size_t e = 0; e < mesh.elements.size(); ++e) {

        Q4NodalValues X{};
        Q4NodalValues u{};

        for (std::size_t i = 0; i < 4; ++i) {

            const std::size_t node =
                mesh.elements[e][i];

            X[i] = mesh.nodes[node];

            u[i] = {
                result.displacement[2 * node],
                result.displacement[2 * node + 1]
            };
        }

        // 7.1 Gauss point checks

        for (double xi : gaussPoints) {

            for (double eta : gaussPoints) {

                const ShapeData shape =
                    computeShapeData(X, xi, eta);

                minDetJ = std::min(minDetJ, shape.detJ);

                const Matrix F =
                    computeIndependentF(u, shape);

                for (std::size_t i = 0; i < 2; ++i) {

                    for (std::size_t j = 0; j < 2; ++j) {

                        maxFError = std::max(
                            maxFError,
                            relativeError(F[i][j], exactF[i][j])
                        );
                    }
                }

                const GreenLagrangeStrain strain =
                    computeGreenLagrangeStrain(F);

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

                    relativeError(stress.sxx, S[0]),
                    relativeError(stress.syy, S[1]),
                    relativeError(stress.sxy, S[2])

                });

                ++gaussCount;
            }
        }

        // 7.2 Production element internal force

        const Q4InternalForce numericalForce =
            computeQ4InternalForce(
                X, u, material, thickness
            );

        // 7.3 Independent analytical integration

        const Q4InternalForce analyticalForce =
            independentInternalForce(X, P);

        double forceScale = 1.0;

        for (std::size_t i = 0; i < 8; ++i) {

            forceScale = std::max(
                forceScale,
                std::abs(analyticalForce[i])
            );
        }

        for (std::size_t i = 0; i < 8; ++i) {

            maxForceError = std::max(
                maxForceError,
                std::abs(
                    numericalForce[i] - analyticalForce[i]
                ) / forceScale
            );
        }

        // 7.4 Global assembly

        for (std::size_t i = 0; i < 4; ++i) {

            const std::size_t node =
                mesh.elements[e][i];

            globalInternalForce[2 * node] +=
                numericalForce[2 * i];

            globalInternalForce[2 * node + 1] +=
                numericalForce[2 * i + 1];
        }

        std::cout << "Element " << e
                  << " verified.\n";
    }

    // ==================================
    // 8. Center equilibrium
    // ==================================

    const double centerResidual = std::hypot(
        globalInternalForce[8],
        globalInternalForce[9]
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

    // ==================================
    // 9. Output
    // ==================================

    std::cout
        << "\n========== Verification ==========\n";

    std::cout << "Final load factor: "
              << result.finalLoadFactor << '\n';

    std::cout << "Gauss points checked: "
              << gaussCount << '\n';

    std::cout << "Minimum det(J): "
              << minDetJ << '\n';

    std::cout << "Maximum F relative error: "
              << maxFError << '\n';

    std::cout << "Maximum strain relative error: "
              << maxStrainError << '\n';

    std::cout << "Maximum PK2 stress relative error: "
              << maxStressError << '\n';

    std::cout << "Maximum element force relative error: "
              << maxForceError << '\n';

    std::cout << "Boundary displacement error: "
              << boundaryError << '\n';

    std::cout << "Center displacement error: "
              << centerError << '\n';

    std::cout << "Center absolute residual: "
              << centerResidual << '\n';

    std::cout << "Center relative residual: "
              << relativeCenterResidual << '\n';

    // ==================================
    // 10. Final checks
    // ==================================

    const bool passed =

        result.converged &&

        std::abs(result.finalLoadFactor - 1.0) < 1.0e-12 &&

        minDetJ > 0.0 &&

        gaussCount == 16 &&

        maxFError < tolerance &&

        maxStrainError < tolerance &&

        maxStressError < tolerance &&

        maxForceError < tolerance &&

        boundaryError < 1.0e-10 &&

        centerError < 1.0e-9 &&

        relativeCenterResidual < 1.0e-8;

    std::cout
        << "\n====================================\n"
        << (passed ? "TEST PASSED" : "TEST FAILED")
        << "\n====================================\n";

    return passed ? 0 : 1;
}
