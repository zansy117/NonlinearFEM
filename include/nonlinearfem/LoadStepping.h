#pragma once

#include "nonlinearfem/Solver.h"

#include <cstddef>
#include <string>
#include <vector>

namespace NonlinearFEM {
    struct LoadSteppingOptions {
        double initialIncrement = 0.1;

        double minimumIncrement = 1.0e-4;

        double maximumIncrement = 0.25;

        double cutbackFactor = 0.5;

        double growthFactor = 1.5;

        std::size_t maxRetries = 8;

        std::size_t maxLoadSteps = 100;

        NewtonRaphsonOptions newtonOptions{};
    };

    struct LoadStepRecord {
        double loadFactor = 0.0;

        double loadIncrement = 0.0;     //增量

        std::size_t iterations = 0;

        std::size_t retries = 0;

        double residualNorm = 0.0;

        bool converged = false;

        std::string failureReason;
    };

    struct LoadSteppingResult {
        GlobalDisplacement displacement;

        double finalLoadFactor = 0.0;

        bool converged = false;

        std::vector<LoadStepRecord> history;

        std::string failureReason;
    };

    [[nodiscard]]
    LoadSteppingResult solveLoadStepping (
        const Q4Mesh& mesh,
        const MaterialParameters& material,
        const NodalLoads& loads,
        const DisplacementConstraints& constraints,
        const GlobalDisplacement& initialdisplacement,
        const LoadSteppingOptions& options,
        double thickness = 1.0
    );
    
}