#include "nonlinearfem/Solver.h"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>

using namespace NonlinearFEM;

int main() {
    Q4Mesh mesh{};
    for (int j = 0; j <= 2; ++j)
        for (int i = 0; i <= 2; ++i)
            mesh.nodes.push_back({0.5 * i, 0.5 * j});
    mesh.elements = {{0,1,4,3}, {1,2,5,4}, {3,4,7,6}, {4,5,8,7}};

    MaterialParameters material{};
    material.youngModulus = 2.1e5;
    material.poissonRatio = 0.3;
    const double thickness = 0.2;
    const double referenceArea = thickness; // Unit reference height.
    const DisplacementConstraints constraints = {
        {0,0,0.0}, {3,0,0.0}, {6,0,0.0}, {0,1,0.0}
    };
    NewtonRaphsonOptions options{};
    options.maxIterations = 30;
    options.residualTolerance = 1.0e-10;
    options.displacementTolerance = 1.0e-12;
    GlobalDisplacement displacement(18, 0.0);
    double previousEndDisplacement = 0.0;

    // Homogeneous uniaxial tension with traction-free transverse direction:
    // Eyy = -nu Exx; lambdaY = sqrt(1 - nu*(lambdaX^2 - 1)).
    // Sxx = Young*Exx; nominal traction Pxx = lambdaX*Sxx.
    // Loads act on the reference edge (dead loads), not the deformed area.
    std::cout << std::scientific << std::setprecision(10);
    std::cout << "stretch,total_force,exact_ux,numerical_ux,max_field_error,free_residual,balance_error\n";
    for (int step = 1; step <= 10; ++step) {
        const double stretch = 1.0 + 0.05 * step;
        const double transverseStretch = std::sqrt(
            1.0 - material.poissonRatio * (stretch * stretch - 1.0));
        const double nominalTraction = material.youngModulus * stretch
                                    * (stretch * stretch - 1.0) / 2.0;
        const double totalForce = nominalTraction * referenceArea;
        const auto externalForce = assembleExternalForce(9, {
            {2,0,totalForce/4.0}, {5,0,totalForce/2.0}, {8,0,totalForce/4.0}
        });
        const auto result = solveNewtonRaphson(
            mesh, material, externalForce, constraints, displacement, options, thickness);
        if (!result.converged) {
            std::cerr << "Step " << step << " failed: " << result.failureReason << '\n';
            return 1;
        }
        displacement = result.displacement;
        const auto system = assembleGlobalSystem(mesh, displacement, material, thickness);
        double maxError = 0.0;
        double residualSquared = 0.0;
        double reactionX = 0.0;
        double reactionY = 0.0;
        for (std::size_t node = 0; node < mesh.nodes.size(); ++node) {
            const auto& X = mesh.nodes[node];
            const double exact[2] = {(stretch-1.0)*X[0], (transverseStretch-1.0)*X[1]};
            for (std::size_t d = 0; d < 2; ++d) {
                const std::size_t dof = 2*node+d;
                if (!std::isfinite(displacement[dof])) return 1;
                maxError = std::max(maxError, std::abs(displacement[dof]-exact[d]));
                const double residual = system.internalForce[dof]-externalForce[dof];
                const bool constrained = (d == 0 && node%3 == 0) || (d == 1 && node == 0);
                if (!constrained) residualSquared += residual*residual;
                else if (d == 0) reactionX += residual;
                else reactionY += residual;
            }
        }
        const double freeResidual = std::sqrt(residualSquared)/totalForce;
        const double balanceError = std::hypot(reactionX+totalForce, reactionY)/totalForce;
        const double endDisplacement = displacement[10];
        std::cout << stretch << ',' << totalForce << ',' << stretch-1.0 << ','
                  << endDisplacement << ',' << maxError << ','
                  << freeResidual << ',' << balanceError << '\n';
        if (!(maxError < 1.0e-8 && freeResidual < 1.0e-8 && balanceError < 1.0e-8
              && endDisplacement > previousEndDisplacement)) {
            std::cerr << "Large-deformation tension verification failed at step " << step << '\n';
            return 1;
        }
        previousEndDisplacement = endDisplacement;
    }
    std::cout << "Large-deformation tension curve verified: 10 load points, up to 50% extension.\n";
    return 0;
}