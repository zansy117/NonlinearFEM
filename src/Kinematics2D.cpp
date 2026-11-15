#include "nonlinearfem/Kinematics2D.h"

namespace NonlinearFEM {

    Matrix2D  computeDeformationGradient (
        const Matrix2D& displacementGradient
    ) {
        Matrix2D deformationGradient{};
        
        for(std::size_t i = 0; i < kSpatialDimension ; ++i ) {
            for(std::size_t j = 0; j < kSpatialDimension ; ++j) {

                deformationGradient[i][j] = displacementGradient[i][j];

                if(i == j) deformationGradient[i][i] += 1.0 ;

            }
        }

        return deformationGradient;
        
    }

    GreenLagrangeStrain computeGreenLagrangeStrain(
        const Matrix2D& deformationGradient
    ) {
        GreenLagrangeStrain Gstrain;

        Gstrain.exx = (deformationGradient[0][0]* deformationGradient[0][0] + deformationGradient[1][0] *deformationGradient[1][0] -1.0) / 2.0;
        Gstrain.eyy = (deformationGradient[1][1]* deformationGradient[1][1] + deformationGradient[0][1] *deformationGradient[0][1] -1.0) / 2.0;
        Gstrain.exy = deformationGradient[0][0]* deformationGradient[0][1] + deformationGradient[1][1] *deformationGradient[1][0];

        return Gstrain ;
    }

}
