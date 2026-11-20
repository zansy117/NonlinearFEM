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

    SecondPKStress_3D computeSecondPKStress_3D (
        const GreenLagrangeStrain_3D& strain ,
        const MaterialParameters& material 
    ) {
        validateMaterialParameters(material);

        if (!std::isfinite(strain.exx) || !std::isfinite(strain.eyy) || !std::isfinite(strain.ezz) || 
            !std::isfinite(strain.exy) || !std::isfinite(strain.eyz) || !std::isfinite(strain.exz)) {
            throw std::invalid_argument(
                "Strain components must be finite."
            );
        }

        SecondPKStress_3D Sstress;

        const double E = material.youngModulus;
        const double nu = material.poissonRatio;
        const double tr = strain.exx + strain.eyy + strain.ezz;
        const double lambda = E * nu /((1.0 + nu) * (1.0 - 2 * nu));
        const double mu = 0.5 * E / (1.0 + nu);

        Sstress.sxx = lambda * tr + 2.0 * mu * strain.exx;
        Sstress.syy = lambda * tr + 2.0 * mu * strain.eyy;
        Sstress.szz = lambda * tr + 2.0 * mu * strain.ezz;

        Sstress.sxy = mu * strain.exy;
        Sstress.sxz = mu * strain.exz;
        Sstress.syz = mu * strain.eyz;

        return Sstress;
    }

    MaterialTangentMatrix_3D computeMaterialTangent_3D (
        const MaterialParameters& material
    ) {
        validateMaterialParameters(material); 

        MaterialTangentMatrix_3D C{};
        
        const double E = material.youngModulus;
        const double nu = material.poissonRatio;
        const double lambda = E * nu /((1.0 + nu) * (1.0 - 2 * nu));
        const double mu = 0.5 * E / (1.0 + nu);

        C[0][0] = C[1][1] = C[2][2] = lambda + 2.0 * mu;
        C[0][1] = C[1][0] = C[0][2] = C[2][0] = C[1][2] = C[2][1] = lambda;
        C[3][3] = C[4][4] = C[5][5] = mu;

        return C;
    }

}