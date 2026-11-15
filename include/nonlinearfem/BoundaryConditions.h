#pragma once

#include "nonlinearfem/Assembly.h"

#include <cstddef>
#include <vector>

namespace NonlinearFEM {
    struct PrescribedDisplacement {
        std::size_t node;
        std::size_t direction;
        double value;
    };

    struct NodalForce {
        std::size_t node;
        std::size_t direction; 
        double value;
    };

    using DisplacementConstraints = std::vector<PrescribedDisplacement>;
    using NodalLoads = std::vector<NodalForce>;

    [[nodiscard]]
    GlobalInternalForce assembleExternalForce(
        std::size_t nodeCount,
        const NodalLoads& loads
    );

    void applyDisplacementConstraints(
        GlobalTangentStiffness& stiffness,
        GlobalInternalForce& residual,
        const GlobalDisplacement& currentDisplacement,
        const DisplacementConstraints& constraints
    );
}