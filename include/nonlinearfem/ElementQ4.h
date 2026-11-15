#pragma once

#include "nonlinearfem/Material.h"

#include <array>
#include <cstddef>

namespace NonlinearFEM {
    constexpr std::size_t kQ4NodeCount = 4;
    constexpr std::size_t kQ4DofCount = 8;

    using Q4NodalValues = 
        std::array<std::array<double , 2> , kQ4NodeCount >; 
     
    using Q4InternalForce = 
        std::array<double , kQ4DofCount>;
          
    [[nodiscard]]
    Q4InternalForce computeQ4InternalForce(
        const Q4NodalValues& referenceCoordinates,
        const Q4NodalValues& nodalDisplacements,
        const MaterialParameters& material,
        double thickness = 1.0
    );

    using Q4TangentStiffness = std::array<std::array<double, kQ4DofCount>,kQ4DofCount>;

    [[nodiscard]]
    Q4TangentStiffness computeQ4TangentStiffness(
        const Q4NodalValues& referenceCoordinates,
        const Q4NodalValues& nodalDisplacements,
        const MaterialParameters& material,
        double thickness = 1.0
    );
    
    using Q4LumpedMass = std::array<double, kQ4DofCount>;

    [[nodiscard]] 
    Q4LumpedMass computeQ4LumpedMass (
        const Q4NodalValues& referenceCoordinates,
        double referenceDensity, //参考密度
        double thickness = 1.0
    );

    [[nodiscard]] 
    double computeQ4StrainEnergy (
        const Q4NodalValues& referenceCoordinates,
        const Q4NodalValues& nodalDisplacements,
        const MaterialParameters& material,
        double thickness = 1.0
    );
}