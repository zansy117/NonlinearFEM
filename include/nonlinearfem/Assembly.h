#pragma once

#include "nonlinearfem/ElementQ4.h"

#include <array>
#include <cstddef>
#include <vector>

namespace NonlinearFEM {
    using Q4NodeCoordinate = std::array<double , 2>;
    
    using Q4ElementConnectivity = std::array<std::size_t , kQ4NodeCount>;
    
    struct Q4Mesh {
        std::vector<Q4NodeCoordinate> nodes;
        std::vector<Q4ElementConnectivity> elements;
    };

    using GlobalDisplacement = std::vector<double>;

    using GlobalInternalForce = std::vector<double>;

    using GlobalTangentStiffness = std::vector<std::vector<double>>;

    struct GlobalAssemblyResult {
        GlobalInternalForce internalForce;
        GlobalTangentStiffness tangentStiffness;
    };

    [[nodiscard]] 
    GlobalAssemblyResult assembleGlobalSystem (
        const Q4Mesh& mesh,
        const GlobalDisplacement& displacement,
        const MaterialParameters& material,
        double thickness = 1.0
    );
    
}