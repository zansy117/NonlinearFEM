#include "nonlinearfem/Material.h" 
#include <cmath>
#include <stdexcept>

namespace {

    void validateMaterialParameters(
        const NonlinearFEM::MaterialParameters& material
    ) {
        if (!std::isfinite(material.youngModulus) ||
            material.youngModulus <= 0.0) {
            throw std::invalid_argument(
                "Young's modulus must be finite and positive."
            );
        }

        if (!std::isfinite(material.poissonRatio) ||
            material.poissonRatio <= -1.0 ||
            material.poissonRatio >= 0.5) {
            throw std::invalid_argument(
                "Poisson ratio must be finite and in (-1.0, 0.5)."
            );
        }
    }

}

namespace NonlinearFEM {
    
    SecondPKStress computeSecondPKStress (
        const GreenLagrangeStrain& strain,
        const MaterialParameters& material 
    ) {
        validateMaterialParameters(material);

        if (!std::isfinite(strain.exx) ||
            !std::isfinite(strain.eyy) ||
            !std::isfinite(strain.exy)) {
            throw std::invalid_argument(
                "Strain components must be finite."
            );
        }

        SecondPKStress Sstress;

        const double E = material.youngModulus;
        const double nu = material.poissonRatio;

        Sstress.sxx = E * (strain.exx + strain.eyy * nu)/(1.0 - nu * nu);
        Sstress.syy = E * (strain.eyy + strain.exx * nu)/(1.0 - nu * nu);
        Sstress.sxy = E * strain.exy / (2.0 * (1.0 + nu));

        return Sstress;
    }

    MaterialTangentMatrix computeMaterialTangent(
        const MaterialParameters& material
    ) {
        validateMaterialParameters(material); 

        MaterialTangentMatrix C{};
        
        const double E = material.youngModulus;
        const double nu = material.poissonRatio;

        C[0][0] = C[1][1] = E / (1.0 - nu * nu);
        C[0][1] = C[1][0] = E * nu / (1.0 - nu * nu);
        C[2][2] = E / (2.0 * (1.0 + nu));

        return C;
    }

}