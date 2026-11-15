
#include "nonlinearfem/Assembly.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <stdexcept>

using namespace NonlinearFEM;

int main() {

    // 1. Construct mesh

    Q4Mesh mesh{};

    mesh.nodes = {
        {{0.0, 0.0}},
        {{1.0, 0.0}},
        {{2.0, 0.0}},
        {{0.0, 1.0}},
        {{1.0, 1.0}},
        {{2.0, 1.0}}
    };

    mesh.elements = {
        {0, 1, 4, 3},
        {1, 2, 5, 4}
    };

    // 2. Material parameters

    const MaterialParameters material{
        2.1e5,
        0.30
    };

    constexpr double thickness = 1.0;

    // 3. Global displacement

    GlobalDisplacement displacement = {
        0.000,  0.000,
        0.010,  0.002,
        0.020,  0.004,
        0.001,  0.005,
        0.012,  0.008,
        0.023,  0.011
    };

    // 4. Assemble global system

    const GlobalAssemblyResult result =
        assembleGlobalSystem(
            mesh,
            displacement,
            material,
            thickness
        );

    // 5. Independently assemble element contributions

    GlobalInternalForce expectedForce(12, 0.0);

    GlobalTangentStiffness expectedStiffness(
        12,
        std::vector<double>(12, 0.0)
    );

    for (const auto& element : mesh.elements) {

        Q4NodalValues coordinates{};
        Q4NodalValues elementDisplacements{};

        std::array<std::size_t, kQ4DofCount> position{};

        for (std::size_t a = 0; a < kQ4NodeCount; ++a) {

            const std::size_t node = element[a];

            coordinates[a] = mesh.nodes[node];

            elementDisplacements[a][0] =
                displacement[2 * node];

            elementDisplacements[a][1] =
                displacement[2 * node + 1];

            position[2 * a] = 2 * node;
            position[2 * a + 1] = 2 * node + 1;
        }

        const auto fe = computeQ4InternalForce(
            coordinates,
            elementDisplacements,
            material,
            thickness
        );

        const auto ke = computeQ4TangentStiffness(
            coordinates,
            elementDisplacements,
            material,
            thickness
        );

        for (std::size_t a = 0; a < kQ4DofCount; ++a) {

            expectedForce[position[a]] += fe[a];

            for (std::size_t b = 0; b < kQ4DofCount; ++b) {

                expectedStiffness[position[a]][position[b]]
                    += ke[a][b];
            }
        }
    }

    // 6. Compare assembled results

    double forceError = 0.0;
    double stiffnessError = 0.0;
    double symmetryError = 0.0;

    for (std::size_t i = 0; i < 12; ++i) {

        forceError = std::max(
            forceError,
            std::abs(result.internalForce[i] - expectedForce[i])
        );

        for (std::size_t j = 0; j < 12; ++j) {

            stiffnessError = std::max(
                stiffnessError,
                std::abs(
                    result.tangentStiffness[i][j]
                    - expectedStiffness[i][j]
                )
            );

            symmetryError = std::max(
                symmetryError,
                std::abs(
                    result.tangentStiffness[i][j]
                    - result.tangentStiffness[j][i]
                )
            );
        }
    }

    // 7. Output

    std::cout << std::scientific << std::setprecision(8);

    std::cout << "Q4 Global Assembly Test\n";
    std::cout << "------------------------------\n";

    std::cout << "Force error     : "
              << forceError << '\n';

    std::cout << "Stiffness error  : "
              << stiffnessError << '\n';

    std::cout << "Symmetry error   : "
              << symmetryError << '\n';

    constexpr double tolerance = 1.0e-10;

    if (forceError < tolerance &&
        stiffnessError < tolerance &&
        symmetryError < tolerance) {

        std::cout << "TEST PASSED\n";
        return 0;
    }

    std::cout << "TEST FAILED\n";
    return 1;
}
