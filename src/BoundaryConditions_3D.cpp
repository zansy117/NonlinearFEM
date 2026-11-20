#include "nonlinearfem/BoundaryConditions_3D.h"

#include <cstddef>
#include <vector>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace NonlinearFEM {
    [[nodiscard]]
    GlobalInternalForce_3D assembleExternalForce_3D (
        std::size_t nodeCount,
        const NodalLoads_3D& loads
    ) {
        if (nodeCount == 0) {
            throw std::invalid_argument(
                "Node count must be positive."
            );
        }

        if (nodeCount > std::numeric_limits<std::size_t>::max() / 3) {
            throw std::overflow_error(
                "Global DOF count overflow."
            );
        }

        for(std::size_t i = 0;i < loads.size();++i) {
            const NodalForce& force = loads[i];

            if (force.node >= nodeCount || force.direction >= 3) {
                throw std::out_of_range(
                    "Invalid nodal load position."
                );
            }

            if (!std::isfinite(force.value)) {
                throw std::invalid_argument(
                    "Nodal load must be finite."
                );
            }
        }

        GlobalInternalForce_3D externalForce;
        externalForce.assign(3 * nodeCount, 0.0);

        for(std::size_t i = 0;i < loads.size(); ++i) {
            NodalForce force = loads[i];
            
            externalForce[3 * force.node + force.direction] += force.value;
        }

        return externalForce;
    }

    void applyDisplacementConstraints_3D(
        GlobalTangentStiffness_3D& stiffness,
        GlobalInternalForce_3D& residual,       
        const GlobalDisplacement_3D& currentDisplacement,
        const DisplacementConstraints_3D& constraints      
    ) {
        const size_t count = currentDisplacement.size();

        if (count == 0 || count % 3 != 0 || stiffness.size() != count || residual.size() != count) {
            throw std::invalid_argument(
                "Invalid three-dimensional system dimensions."
            );
        }

        for(std::size_t i = 0;i < count;++i) {
            if (stiffness[i].size() != count) {
                throw std::invalid_argument(
                    "Global stiffness matrix must be square."
                );
            }

            if (!std::isfinite(currentDisplacement[i]) || !std::isfinite(residual[i])) {
                throw std::invalid_argument(
                    "Displacement and residual must be finite."
                );
            }

            for(std::size_t j = 0;j < count;++j) {
                if (!std::isfinite(stiffness[i][j])) {
                    throw std::invalid_argument(
                        "Global stiffness entries must be finite."
                    );
                }
            }
        }

        std::vector<bool> constrained(count,false);

        for(std::size_t i = 0;i < constraints.size();++i) {
            const PrescribedDisplacement& constraint = constraints[i];

            if (constraint.node >= count / 3 || constraint.direction >= 3) {
                throw std::out_of_range(
                    "Invalid displacement constraint position."
                );
            }

            if (!std::isfinite(constraint.value)) {
                throw std::invalid_argument(
                    "Prescribed displacement must be finite."
                );
            }

            const size_t position = 3 * constraint.node + constraint.direction;

            if (constrained[position]) {
                throw std::invalid_argument(
                    "Duplicate displacement constraint."
                );
            }

            if (!std::isfinite(constraint.value - currentDisplacement[position])) {
                throw std::overflow_error(
                    "Displacement increment overflow."
                );
            }

            constrained[position] = true;
        }

        for(std::size_t i = 0;i < constraints.size() ;++i) {
            const size_t position = constraints[i].node * 3 + constraints[i].direction;

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