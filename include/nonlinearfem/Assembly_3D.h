#pragma once

#include "nonlinearfem/ElementH8.h"

#include <array>
#include <cstddef>
#include <vector>

namespace NonlinearFEM {
    using H8NodeCoordinate = std::array<double , 3>;
    
    using H8ElementConnectivity = std::array<std::size_t , kH8NodeCount>;
    
    struct H8Mesh {
        std::vector<H8NodeCoordinate> nodes;
        std::vector<H8ElementConnectivity> elements;
    };

    using GlobalDisplacement_3D = std::vector<double>;

    using GlobalInternalForce_3D = std::vector<double>;

    using GlobalTangentStiffness_3D = std::vector<std::vector<double>>;

    struct GlobalAssemblyResult_3D {
        GlobalInternalForce_3D internalForce;
        GlobalTangentStiffness_3D tangentStiffness;
    };

    [[nodiscard]]
    GlobalAssemblyResult_3D assembleGlobalSystem_3D (
        const H8Mesh& mesh,
        const GlobalDisplacement_3D& displacement,
        const MaterialParameters& material
    );

    [[nodiscard]]
    GlobalInternalForce_3D assembleInternalForce_3D (
        const H8Mesh& mesh,
        const GlobalDisplacement_3D& displacement,
        const MaterialParameters& material
    );
    
    using GlobalLumpedMass_3D = std::vector<double>;

    [[nodiscard]]
    GlobalLumpedMass_3D assembleGlobalLumpedMass_3D (
        const H8Mesh& mesh,
        double referenceDensity
    );

    double assembleGlobalEnergy_3D (
        const H8Mesh& mesh,
        const GlobalDisplacement_3D& displacement,
        const MaterialParameters& material
    );

}