
#include "nonlinearfem/Assembly.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace NonlinearFEM;

namespace {

constexpr double kThickness = 1.0;
constexpr double kTolerance = 1.0e-10;

int passedTests = 0;
int failedTests = 0;


// ============================================================
// Test utilities
// ============================================================

void reportTest(
    const std::string& name,
    bool passed,
    double error = 0.0
) {
    if (passed) {
        ++passedTests;
        std::cout << "[PASSED] " << name
                  << " | Error: " << error << '\n';
    }
    else {
        ++failedTests;
        std::cout << "[FAILED] " << name
                  << " | Error: " << error << '\n';
    }
}


double computeVectorNorm(
    const std::vector<double>& v
) {
    double sum = 0.0;

    for (double value : v) {
        sum += value * value;
    }

    return std::sqrt(sum);
}


double computeMatrixNorm(
    const std::vector<std::vector<double>>& A
) {
    double sum = 0.0;

    for (const auto& row : A) {
        for (double value : row) {
            sum += value * value;
        }
    }

    return std::sqrt(sum);
}


double relativeVectorError(
    const std::vector<double>& a,
    const std::vector<double>& b
) {
    if (a.size() != b.size()) {
        throw std::invalid_argument(
            "Vector dimensions do not match."
        );
    }

    std::vector<double> difference(a.size());

    for (std::size_t i = 0; i < a.size(); ++i) {
        difference[i] = a[i] - b[i];
    }

    return computeVectorNorm(difference) /
        std::max(computeVectorNorm(b), 1.0e-30);
}


double relativeMatrixError(
    const std::vector<std::vector<double>>& A,
    const std::vector<std::vector<double>>& B
) {
    if (A.size() != B.size()) {
        throw std::invalid_argument(
            "Matrix dimensions do not match."
        );
    }

    double differenceSquared = 0.0;
    double referenceSquared = 0.0;

    for (std::size_t i = 0; i < A.size(); ++i) {

        if (A[i].size() != B[i].size()) {
            throw std::invalid_argument(
                "Matrix dimensions do not match."
            );
        }

        for (std::size_t j = 0; j < A[i].size(); ++j) {
            const double difference = A[i][j] - B[i][j];

            differenceSquared += difference * difference;
            referenceSquared += B[i][j] * B[i][j];
        }
    }

    return std::sqrt(differenceSquared) /
        std::max(std::sqrt(referenceSquared), 1.0e-30);
}


// ============================================================
// Test mesh and material
// ============================================================

Q4Mesh createTestMesh() {
    Q4Mesh mesh;

    mesh.nodes = {
        {0.0, 0.0},
        {1.0, 0.0},
        {2.1, 0.1},
        {0.0, 1.0},
        {1.1, 1.1},
        {2.2, 1.2}
    };

    mesh.elements = {
        {0, 1, 4, 3},
        {1, 2, 5, 4}
    };

    return mesh;
}


MaterialParameters createTestMaterial() {
    MaterialParameters material;

    material.youngModulus = 2.1e5;
    material.poissonRatio = 0.3;

    return material;
}


GlobalDisplacement createTestDisplacement() {
    return {
         0.0000,  0.0000,
         0.0010,  0.0003,
         0.0018,  0.0008,
        -0.0002,  0.0007,
         0.0008,  0.0015,
         0.0015,  0.0020
    };
}


// ============================================================
// Independent element-to-global scatter
// ============================================================

GlobalAssemblyResult independentlyScatterElements(
    const Q4Mesh& mesh,
    const GlobalDisplacement& displacement,
    const MaterialParameters& material
) {
    const std::size_t dofCount = 2 * mesh.nodes.size();

    GlobalAssemblyResult result{};

    result.internalForce.assign(dofCount, 0.0);

    result.tangentStiffness.assign(
        dofCount,
        std::vector<double>(dofCount, 0.0)
    );

    for (const auto& element : mesh.elements) {

        Q4NodalValues coordinates{};
        Q4NodalValues localDisplacements{};

        std::array<std::size_t, kQ4DofCount> dofMap{};

        for (std::size_t a = 0; a < kQ4NodeCount; ++a) {

            const std::size_t node = element[a];

            coordinates[a] = mesh.nodes[node];

            localDisplacements[a][0] =
                displacement[2 * node];

            localDisplacements[a][1] =
                displacement[2 * node + 1];

            dofMap[2 * a] = 2 * node;
            dofMap[2 * a + 1] = 2 * node + 1;
        }

        const auto localForce = computeQ4InternalForce(
            coordinates,
            localDisplacements,
            material,
            kThickness
        );

        const auto localStiffness = computeQ4TangentStiffness(
            coordinates,
            localDisplacements,
            material,
            kThickness
        );

        // Scatter force
        for (std::size_t a = 0; a < kQ4DofCount; ++a) {
            result.internalForce[dofMap[a]] += localForce[a];
        }

        // Scatter stiffness
        for (std::size_t a = 0; a < kQ4DofCount; ++a) {
            for (std::size_t b = 0; b < kQ4DofCount; ++b) {

                result.tangentStiffness[dofMap[a]][dofMap[b]]
                    += localStiffness[a][b];
            }
        }
    }

    return result;
}


// ============================================================
// Test 5.1 and Test 5.2
// ============================================================

void testGlobalAssembly() {
    const Q4Mesh mesh = createTestMesh();

    const auto displacement = createTestDisplacement();

    const auto material = createTestMaterial();

    const auto assembled = assembleGlobalSystem(
        mesh,
        displacement,
        material,
        kThickness
    );

    const auto reference = independentlyScatterElements(
        mesh,
        displacement,
        material
    );

    const double forceError = relativeVectorError(
        assembled.internalForce,
        reference.internalForce
    );

    const double stiffnessError = relativeMatrixError(
        assembled.tangentStiffness,
        reference.tangentStiffness
    );

    reportTest(
        "Test 5.1 - Global Internal Force Assembly",
        forceError < kTolerance,
        forceError
    );

    reportTest(
        "Test 5.2 - Global Tangent Stiffness Assembly",
        stiffnessError < kTolerance,
        stiffnessError
    );
}


// ============================================================
// Test 5.3
// Global tangent finite difference verification
// ============================================================

void testGlobalTangentFiniteDifference() {
    const Q4Mesh mesh = createTestMesh();

    const auto displacement = createTestDisplacement();

    const auto material = createTestMaterial();

    const auto analytical = assembleGlobalSystem(
        mesh,
        displacement,
        material,
        kThickness
    );

    const std::size_t dofCount = displacement.size();

    constexpr double step = 1.0e-6;

    std::vector<std::vector<double>> numerical(
        dofCount,
        std::vector<double>(dofCount, 0.0)
    );

    // Central finite difference for every global DOF
    for (std::size_t j = 0; j < dofCount; ++j) {

        GlobalDisplacement displacementPlus = displacement;
        GlobalDisplacement displacementMinus = displacement;

        displacementPlus[j] += step;
        displacementMinus[j] -= step;

        const auto forcePlus = assembleGlobalSystem(
            mesh,
            displacementPlus,
            material,
            kThickness
        );

        const auto forceMinus = assembleGlobalSystem(
            mesh,
            displacementMinus,
            material,
            kThickness
        );

        for (std::size_t i = 0; i < dofCount; ++i) {

            numerical[i][j] =
                (forcePlus.internalForce[i] -
                 forceMinus.internalForce[i]) /
                (2.0 * step);
        }
    }

    const double error = relativeMatrixError(
        analytical.tangentStiffness,
        numerical
    );

    reportTest(
        "Test 5.3 - Global Tangent Finite Difference",
        error < 1.0e-6,
        error
    );
}


// ============================================================
// Test 5.4
// Element ordering invariance
// ============================================================

void testElementOrderingInvariance() {
    Q4Mesh meshA = createTestMesh();
    Q4Mesh meshB = meshA;

    std::reverse(
        meshB.elements.begin(),
        meshB.elements.end()
    );

    const auto displacement = createTestDisplacement();

    const auto material = createTestMaterial();

    const auto resultA = assembleGlobalSystem(
        meshA,
        displacement,
        material,
        kThickness
    );

    const auto resultB = assembleGlobalSystem(
        meshB,
        displacement,
        material,
        kThickness
    );

    const double forceError = relativeVectorError(
        resultA.internalForce,
        resultB.internalForce
    );

    const double stiffnessError = relativeMatrixError(
        resultA.tangentStiffness,
        resultB.tangentStiffness
    );

    constexpr double orderingTolerance = 1.0e-12;

    reportTest(
        "Test 5.4 - Element Ordering Invariance",
        forceError < orderingTolerance &&
        stiffnessError < orderingTolerance,
        std::max(forceError, stiffnessError)
    );
}

} // anonymous namespace


// ============================================================
// Main
// ============================================================

int main() {
    try {
        std::cout
            << "========================================\n"
            << "Multi-element Global Assembly Test\n"
            << "========================================\n\n";

        std::cout << "Nodes: 6\n";
        std::cout << "Elements: 2\n";
        std::cout << "DOFs: 12\n\n";

        testGlobalAssembly();

        testGlobalTangentFiniteDifference();

        testElementOrderingInvariance();
    }
    catch (const std::exception& error) {

        std::cerr << "\n[FATAL ERROR] "
                  << error.what() << '\n';

        return 2;
    }

    std::cout
        << "\n========================================\n"
        << "Verification Summary\n"
        << "========================================\n";

    std::cout << "Passed: " << passedTests << '\n';
    std::cout << "Failed: " << failedTests << '\n';

    if (failedTests == 0) {
        std::cout
            << "\nALL TESTS PASSED\n";
        return 0;
    }

    std::cout
        << "\nTEST FAILED\n";

    return 1;
}
