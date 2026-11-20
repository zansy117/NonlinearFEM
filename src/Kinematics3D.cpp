#include "nonlinearfem/Kinematics3D.h"

namespace NonlinearFEM {

    Matrix3D  computeDeformationGradient_3D (
        const Matrix3D& displacementGradient
    ) {
        Matrix3D deformationGradient{};
        
        for(std::size_t i = 0; i < kSpatialDimension_3D ; ++i ) {
            for(std::size_t j = 0; j < kSpatialDimension_3D ; ++j) {
                deformationGradient[i][j] = displacementGradient[i][j];

                if(i == j) deformationGradient[i][i] += 1.0 ;

            }
        }

        return deformationGradient;
        
    }

    GreenLagrangeStrain_3D computeGreenLagrangeStrain_3D(
        const Matrix3D& deformationGradient
    ) {
        GreenLagrangeStrain_3D Gstrain;

        Gstrain.exx = (deformationGradient[0][0]* deformationGradient[0][0] + deformationGradient[1][0] *deformationGradient[1][0] + deformationGradient[2][0] * deformationGradient[2][0] - 1.0) / 2.0;
        Gstrain.eyy = (deformationGradient[0][1]* deformationGradient[0][1] + deformationGradient[1][1] *deformationGradient[1][1] + deformationGradient[2][1] * deformationGradient[2][1] - 1.0) / 2.0;
        Gstrain.ezz = (deformationGradient[0][2]* deformationGradient[0][2] + deformationGradient[1][2] *deformationGradient[1][2] + deformationGradient[2][2] * deformationGradient[2][2] - 1.0) / 2.0;
        Gstrain.exy = deformationGradient[0][0]* deformationGradient[0][1] + deformationGradient[1][0] *deformationGradient[1][1] + deformationGradient[2][0] * deformationGradient[2][1];
        Gstrain.exz = deformationGradient[0][0]* deformationGradient[0][2] + deformationGradient[1][0] *deformationGradient[1][2] + deformationGradient[2][0] * deformationGradient[2][2];
        Gstrain.eyz = deformationGradient[0][1]* deformationGradient[0][2] + deformationGradient[1][1] *deformationGradient[1][2] + deformationGradient[2][1] * deformationGradient[2][2];

        return Gstrain ;
    }

}
