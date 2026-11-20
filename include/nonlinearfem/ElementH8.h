#pragma once

#include "nonlinearfem/Material.h"

#include <array>
#include <cstddef>

namespace NonlinearFEM {
    constexpr std::size_t kH8NodeCount = 8;
    constexpr std::size_t kH8DofCount = 24;

    using H8NodalValues = 
        std::array<std::array<double , 3> , kH8NodeCount >; 
     
    using H8InternalForce = 
        std::array<double , kH8DofCount>;
          
    [[nodiscard]]
    H8InternalForce computeH8InternalForce(
        const H8NodalValues& referenceCoordinates,
        const H8NodalValues& nodalDisplacements,
        const MaterialParameters& material
    );

    using H8TangentStiffness = std::array<std::array<double, kH8DofCount>,kH8DofCount>;

    [[nodiscard]]
    H8TangentStiffness computeH8TangentStiffness(
        const H8NodalValues& referenceCoordinates,
        const H8NodalValues& nodalDisplacements,
        const MaterialParameters& material
    );
    
    using H8LumpedMass = std::array<double, kH8DofCount>;

    [[nodiscard]] 
    H8LumpedMass computeH8LumpedMass (
        const H8NodalValues& referenceCoordinates,
        double referenceDensity
    );

    [[nodiscard]] 
    double computeH8StrainEnergy (
        const H8NodalValues& referenceCoordinates,
        const H8NodalValues& nodalDisplacements,
        const MaterialParameters& material
    );
}