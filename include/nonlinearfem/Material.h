#pragma once

#include "nonlinearfem/Kinematics2D.h"
#include <array>

namespace NonlinearFEM {

    struct MaterialParameters {
        double youngModulus = 0.0;
        double poissonRatio = 0.0;
    }; 

    struct SecondPKStress {
        double sxx = 0.0;
        double syy = 0.0;
        double sxy = 0.0;
    };

    constexpr std::size_t kVoigtDimension = 3;

    using MaterialTangentMatrix =
        std::array<
            std::array<double, kVoigtDimension>,
            kVoigtDimension
        >;
    
    [[nodiscard]] 
    SecondPKStress computeSecondPKStress (
        const GreenLagrangeStrain& strain ,
        const MaterialParameters& material 
    );

    [[nodiscard]] 
    MaterialTangentMatrix computeMaterialTangent(
        const MaterialParameters& material
    );

} // namespace NonlinearFEM