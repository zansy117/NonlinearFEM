#pragma once

#include "nonlinearfem/Assembly_3D.h"
#include "nonlinearfem/BoundaryConditions.h"


namespace NonlinearFEM {
    using DisplacementConstraints_3D = std::vector<PrescribedDisplacement>;
    using NodalLoads_3D = std::vector<NodalForce>;

    [[nodiscard]]
    GlobalInternalForce_3D assembleExternalForce_3D(
        std::size_t nodeCount,
        const NodalLoads_3D& loads
    );

    void applyDisplacementConstraints_3D(
        GlobalTangentStiffness_3D& stiffness,
        GlobalInternalForce_3D& residual,
        const GlobalDisplacement_3D& currentDisplacement,
        const DisplacementConstraints_3D& constraints
    );
}