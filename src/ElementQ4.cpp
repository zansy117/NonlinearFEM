#include "nonlinearfem/ElementQ4.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>


namespace NonlinearFEM {

    namespace {
        struct Q4ShapeFunctions {
            std::array<double,4> N{};
            std::array<double,4> dNdxi{};
            std::array<double,4> dNdeta{};
        };

        [[nodiscard]]
        Q4ShapeFunctions computeQ4Shapefunctions (double xi,double eta);
        Q4ShapeFunctions computeQ4Shapefunctions (double xi,double eta) {
            Q4ShapeFunctions shapefunction ;
            shapefunction.N[0]=(1.0 - xi) * (1.0 - eta) / 4.0;
            shapefunction.N[1]=(1.0 + xi) * (1.0 - eta) / 4.0;
            shapefunction.N[2]=(1.0 + xi) * (1.0 + eta) / 4.0;
            shapefunction.N[3]=(1.0 - xi) * (1.0 + eta) / 4.0;
        
            shapefunction.dNdxi[0]=(eta - 1.0) / 4.0;
            shapefunction.dNdxi[1]=(1.0 - eta) / 4.0;
            shapefunction.dNdxi[2]=(eta + 1.0) / 4.0;
            shapefunction.dNdxi[3]=(-eta - 1.0) / 4.0;

            shapefunction.dNdeta[0]=(xi - 1.0) / 4.0;
            shapefunction.dNdeta[1]=(-1.0 - xi) / 4.0;
            shapefunction.dNdeta[2]=(xi + 1.0) / 4.0;
            shapefunction.dNdeta[3]=(1.0 - xi) / 4.0;

            return shapefunction ;
        }

        struct JacobianData {
            std::array<std::array<double, 2>, 2> JacobianMatrix{};

            double detJ=0.0;

            std::array<std::array<double, 2>, 4> dNdX{};
        };

        [[nodiscard]]
        JacobianData computeQ4Jacobian(
            const Q4NodalValues& referenceCoordinates,
            const Q4ShapeFunctions& shape
        ) {
            JacobianData J{};

            for(std::size_t i = 0;i < kQ4NodeCount ; ++i) {
                J.JacobianMatrix[0][0] += shape.dNdxi[i] * referenceCoordinates[i][0];
                J.JacobianMatrix[0][1] += shape.dNdeta[i] * referenceCoordinates[i][0];
                J.JacobianMatrix[1][0] += shape.dNdxi[i] * referenceCoordinates[i][1];
                J.JacobianMatrix[1][1] += shape.dNdeta[i] * referenceCoordinates[i][1];
            }

            J.detJ = J.JacobianMatrix[0][0] * J.JacobianMatrix[1][1] - J.JacobianMatrix[1][0] * J.JacobianMatrix[0][1];

            const double scale = std::max(std::max(std::abs(J.JacobianMatrix[0][0]),std::abs(J.JacobianMatrix[0][1])),std::max(std::abs(J.JacobianMatrix[1][0]),std::abs(J.JacobianMatrix[1][1])));

            const double tolerance = 100.0 * std::numeric_limits<double>::epsilon() * scale * scale;

            if (!std::isfinite(J.detJ) || J.detJ <= tolerance) {
                throw std::invalid_argument(
                    "Invalid Q4 Jacobian."
                );
            }

            std::array<std::array<double , 2>, 2> invJ{};
            invJ[0][0] = J.JacobianMatrix[1][1]/J.detJ;
            invJ[1][0] = -J.JacobianMatrix[1][0]/J.detJ;
            invJ[0][1] = -J.JacobianMatrix[0][1]/J.detJ;
            invJ[1][1] = J.JacobianMatrix[0][0]/J.detJ;

            for(std::size_t i = 0;i < kQ4NodeCount ; ++i) {
                J.dNdX[i][0] = shape.dNdxi[i] * invJ[0][0] + shape.dNdeta[i] * invJ[1][0];
                J.dNdX[i][1] = shape.dNdxi[i] * invJ[0][1] + shape.dNdeta[i] * invJ[1][1];
            }

            return J;

        }

        Matrix2D computeQ4DisplacementGradient (
            const Q4NodalValues& nodaldisplacements,
            const JacobianData& jacobian 
        )  {
            Matrix2D displacementgradient{};

            for(std::size_t i = 0 ;i < kQ4NodeCount ; ++i) {
                displacementgradient[0][0] += jacobian.dNdX[i][0] * nodaldisplacements[i][0]; 
                displacementgradient[0][1] += jacobian.dNdX[i][1] * nodaldisplacements[i][0]; 
                displacementgradient[1][0] += jacobian.dNdX[i][0] * nodaldisplacements[i][1]; 
                displacementgradient[1][1] += jacobian.dNdX[i][1] * nodaldisplacements[i][1];
            }

            return displacementgradient;
        }

        struct Q4GaussPointState {
            double detJ = 0.0;
            Matrix2D deformationGradient{};
            GreenLagrangeStrain strain{};
            SecondPKStress stress{};
            std::array<std::array<double, 2>, 4> dNdX{};
        };

        [[nodiscard]]
        Q4GaussPointState computeQ4GaussPointState (
            double xi,
            double eta,
            const Q4NodalValues& referenceCoordinates,
            const Q4NodalValues& nodaldisplacements,
            const MaterialParameters& material
        ) {
            const Q4ShapeFunctions shape = computeQ4Shapefunctions(xi , eta);

            const JacobianData jacobian = computeQ4Jacobian(referenceCoordinates , shape);

            const Matrix2D displacementgradient = computeQ4DisplacementGradient(nodaldisplacements , jacobian);

            const Matrix2D deformationgradient = computeDeformationGradient(displacementgradient);

            const GreenLagrangeStrain strain = computeGreenLagrangeStrain(deformationgradient);

            const SecondPKStress stress = computeSecondPKStress(strain , material);

            Q4GaussPointState state{};
            state.detJ = jacobian.detJ;
            state.deformationGradient = deformationgradient;
            state.strain = strain;
            state.stress = stress;
            state.dNdX = jacobian.dNdX;
            
            return state;
        }

        [[nodiscard]]
        Matrix2D computeFirstPKStress (
            const Matrix2D& deformationgradient,
            const SecondPKStress& stress
        ) {
            Matrix2D firstPKStress{};

            firstPKStress[0][0] = deformationgradient[0][0] * stress.sxx + deformationgradient[0][1] * stress.sxy;
            firstPKStress[0][1] = deformationgradient[0][0] * stress.sxy + deformationgradient[0][1] * stress.syy;
            firstPKStress[1][0] = deformationgradient[1][0] * stress.sxx + deformationgradient[1][1] * stress.sxy;
            firstPKStress[1][1] = deformationgradient[1][0] * stress.sxy + deformationgradient[1][1] * stress.syy;

            return firstPKStress;
        }
    }



    Q4InternalForce computeQ4InternalForce(
        const Q4NodalValues& referenceCoordinates,
        const Q4NodalValues& nodalDisplacements,
        const MaterialParameters& material,
        double thickness
    ) {
        Q4InternalForce internalForce{};

        constexpr double gauss = 0.57735026918962576451;

        const std::array<double , 2> gausspoints = {-gauss,gauss};

        for(std::size_t i = 0 ;i < 2; ++i) {
            for(std::size_t j = 0; j < 2; ++j) {
                const double xi = gausspoints[i] ,eta = gausspoints[j];

                const Q4GaussPointState state = computeQ4GaussPointState(xi , eta , referenceCoordinates , nodalDisplacements ,material);

                const Matrix2D firstPKStress = computeFirstPKStress(state.deformationGradient,state.stress);

                const double integrationWeight = thickness * state.detJ;

                for(std::size_t k = 0; k < kQ4NodeCount ; ++k) {
                    double forcex = 0.0, forcey = 0.0;

                    for(std::size_t h = 0;h < 2 ; ++h) {
                        forcex += state.dNdX[k][h] * firstPKStress[0][h];
                        forcey += state.dNdX[k][h] * firstPKStress[1][h];
                    } 

                    internalForce[2*k] += forcex * integrationWeight;
                    internalForce[2*k+1] += forcey * integrationWeight;
                }
            }
        }

        return internalForce;
    }

    [[nodiscard]]
    Q4TangentStiffness computeQ4TangentStiffness(
        const Q4NodalValues& referenceCoordinates,
        const Q4NodalValues& nodalDisplacements,
        const MaterialParameters& material,
        double thickness 
    ) {
        Q4TangentStiffness K{} ;

        const MaterialTangentMatrix C = computeMaterialTangent(material);

        constexpr double gauss = 0.57735026918962576451;

        const std::array<double , 2> gausspoints = {-gauss,gauss};

        std::array<std::array<double,kQ4DofCount>,3> B{};

        for(std::size_t i = 0; i < 2; ++i) {
            for(std::size_t j = 0 ;j < 2; ++j) {
                const double xi = gausspoints[i],eta = gausspoints[j];

                Q4GaussPointState state = computeQ4GaussPointState(xi,eta,referenceCoordinates,nodalDisplacements,material);

                for(std::size_t k = 0 ;k < kQ4NodeCount ;++k) {
                    B[0][2*k] = state.deformationGradient[0][0] * state.dNdX[k][0];
                    B[0][2*k+1] = state.deformationGradient[1][0] * state.dNdX[k][0];
                    B[1][2*k] = state.deformationGradient[0][1] * state.dNdX[k][1];
                    B[1][2*k+1] = state.deformationGradient[1][1] * state.dNdX[k][1];
                    B[2][2*k] = state.deformationGradient[0][1] * state.dNdX[k][0] + state.deformationGradient[0][0] *state.dNdX[k][1];
                    B[2][2*k+1] = state.deformationGradient[1][1] * state.dNdX[k][0] + state.deformationGradient[1][0] * state.dNdX[k][1];
                }

                const double integrationWeight = thickness * state.detJ;


                for(std::size_t k = 0 ;k < kQ4DofCount ;++k) {
                    for(std::size_t h = 0; h < kQ4DofCount ;++h) {
                        for(size_t r = 0;r < 3 ;++r) {
                            for(size_t s = 0 ;s < 3;++s) {
                                K[k][h] += B[r][k] * C[r][s] * B[s][h] *integrationWeight;
                            }
                        }
                    }
                }

                const std::array<std::array<double,2> , 2> S={{
                    {state.stress.sxx , state.stress.sxy} ,
                    {state.stress.sxy , state.stress.syy} 
                }};

                for(std::size_t k = 0 ;k < kQ4NodeCount ;++k) {
                    for(std::size_t h = 0; h < kQ4NodeCount ;++h) {
                        for(size_t r = 0;r < 2 ;++r) {
                            for(size_t s = 0 ;s < 2;++s) {
                                K[k*2][h*2] += state.dNdX[k][r] * S[r][s] * state.dNdX[h][s] * integrationWeight;
                                K[k*2+1][h*2+1] += state.dNdX[k][r] * S[r][s] * state.dNdX[h][s] * integrationWeight;
                            }
                        }
                    }
                }
            }
        }
        return K;
    }
}