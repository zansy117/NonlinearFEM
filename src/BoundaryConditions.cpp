#include "nonlinearfem/BoundaryConditions.h"

#include <cstddef>
#include <vector>

namespace NonlinearFEM {
    [[nodiscard]]
    GlobalInternalForce assembleExternalForce (
        std::size_t nodeCount,
        const NodalLoads& loads
    ) {
        GlobalInternalForce externalForce;
        externalForce.assign(2 * nodeCount, 0.0);

        for(std::size_t i = 0;i < loads.size(); ++i) {
            NodalForce force = loads[i];
            
            externalForce[2 * force.node + force.direction] += force.value;
        }

        return externalForce;
    }

    void applyDisplacementConstraints(
        GlobalTangentStiffness& stiffness,
        GlobalInternalForce& residual,       //残差
        const GlobalDisplacement& currentDisplacement,
        const DisplacementConstraints& constraints      //限制
    ) {
        const size_t count = currentDisplacement.size();

        std::vector<bool> constrained(count,false);

        for(std::size_t i = 0;i < constraints.size() ;++i) {
            const size_t position = constraints[i].node * 2 + constraints[i].direction;

            const double delta = constraints[i].value - currentDisplacement[position];

            constrained[position] = 1;

            for(size_t j = 0;j < count;++j) {
                if(position != j) residual[j] -= stiffness[j][position] * delta ;
            }

            for(std::size_t j = 0;j < count;++j) {
                stiffness[position][j] = 0.0;
                stiffness[j][position] = 0.0;
            }

            stiffness[position][position] = 1.0;

            residual[position] = delta;
        }

        return ;
    }
}