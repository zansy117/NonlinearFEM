#pragma once

#include"nonlinearfem/BoundaryConditions.h"

#include <cstddef>
#include <vector>
#include<string>

namespace NonlinearFEM {
    struct NewtonRaphsonOptions {
        std::size_t maxIterations = 30;  //最大迭代次数

        double residualTolerance = 1.0e-8;

        double displacementTolerance = 1.0e-10;

    };

    struct NewtonRaphsonResult {
        GlobalDisplacement displacement;

        bool converged = false;  //收敛

        std::size_t iterations = 0;     

        double residualNorm = 0.0;

        double displacementIncrementNorm = 0.0;

        std::string failureReason;
    };

    [[nodiscard]]
    NewtonRaphsonResult solveNewtonRaphson(
        const Q4Mesh& mesh,
        const MaterialParameters& material,
        const GlobalInternalForce& externalForce,
        const DisplacementConstraints& constraints,
        const GlobalDisplacement& initialDisplacement,
        const NewtonRaphsonOptions& options,
        double thickness = 1.0
    );


}