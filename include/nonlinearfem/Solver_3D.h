#pragma once

#include"nonlinearfem/BoundaryConditions_3D.h"
#include"nonlinearfem/Solver.h"

#include <cstddef>
#include <vector>
#include<string>

namespace NonlinearFEM {
    struct NewtonRaphsonResult_3D {
        GlobalDisplacement_3D displacement;

        bool converged = false;  

        std::size_t iterations = 0;     

        double residualNorm = 0.0;

        double displacementIncrementNorm = 0.0;

        std::string failureReason;
    };

    [[nodiscard]]
    NewtonRaphsonResult_3D solveNewtonRaphson_3D(
        const H8Mesh& mesh,
        const MaterialParameters& material,
        const GlobalInternalForce_3D& externalForce,
        const DisplacementConstraints_3D& constraints,
        const GlobalDisplacement_3D& initialDisplacement,
        const NewtonRaphsonOptions& options
    );


}