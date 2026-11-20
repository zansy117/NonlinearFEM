#pragma once

#include "nonlinearfem/Solver_3D.h"
#include "nonlinearfem/LoadStepping.h"

#include <cstddef>
#include <string>
#include <vector>

namespace NonlinearFEM {
    struct LoadSteppingResult_3D {
        GlobalDisplacement_3D displacement;

        double finalLoadFactor = 0.0;

        bool converged = false;

        std::vector<LoadStepRecord> history;

        std::string failureReason;
    };

    [[nodiscard]]
    LoadSteppingResult_3D solveLoadStepping_3D (
        const H8Mesh& mesh,
        const MaterialParameters& material,
        const NodalLoads_3D& loads,
        const DisplacementConstraints_3D& constraints,
        const GlobalDisplacement_3D& initialdisplacement,
        const LoadSteppingOptions& options
    );
    
}