#pragma once

#include "nonlinearfem/Kinematics2D.h"
#include "nonlinearfem/Kinematics3D.h"
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

    struct SecondPKStress_3D {
        double sxx = 0.0;
        double syy = 0.0;
        double szz = 0.0;
        double sxy = 0.0;
        double sxz = 0.0;
        double syz = 0.0;
    };

    constexpr std::size_t kVoigtDimension_3D = 6;

    using MaterialTangentMatrix_3D = std::array<std::array<double, kVoigtDimension_3D>,kVoigtDimension_3D>;

    [[nodiscard]] 
    SecondPKStress_3D computeSecondPKStress_3D (
        const GreenLagrangeStrain_3D& strain ,
        const MaterialParameters& material 
    );

    [[nodiscard]] 
    MaterialTangentMatrix_3D computeMaterialTangent_3D (
        const MaterialParameters& material
    );

} // namespace NonlinearFEM