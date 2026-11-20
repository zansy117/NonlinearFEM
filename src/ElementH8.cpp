#include "nonlinearfem/ElementH8.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace NonlinearFEM {
    namespace {
        struct H8ShapeFunctions {
            std::array<double,8> N{};
            std::array<double,8> dNdxi{};
            std::array<double,8> dNdeta{};
            std::array<double,8> dNdzeta{};
        };

        [[nodiscard]]
        H8ShapeFunctions computeH8Shapefunctions (double xi,double eta,double zeta);
        H8ShapeFunctions computeH8Shapefunctions (double xi,double eta,double zeta) {
            H8ShapeFunctions shapefunction ;

            shapefunction.N[0]=(1.0 - xi) * (1.0 - eta) * (1.0 - zeta) / 8.0;
            shapefunction.N[1]=(1.0 - xi) * (1.0 - eta) * (1.0 + zeta) / 8.0;
            shapefunction.N[2]=(1.0 + xi) * (1.0 - eta) * (1.0 + zeta) / 8.0;
            shapefunction.N[3]=(1.0 + xi) * (1.0 - eta) * (1.0 - zeta) / 8.0;
            shapefunction.N[4]=(1.0 + xi) * (1.0 + eta) * (1.0 - zeta) / 8.0;
            shapefunction.N[5]=(1.0 - xi) * (1.0 + eta) * (1.0 - zeta) / 8.0;
            shapefunction.N[6]=(1.0 - xi) * (1.0 + eta) * (1.0 + zeta) / 8.0;
            shapefunction.N[7]=(1.0 + xi) * (1.0 + eta) * (1.0 + zeta) / 8.0;

        
            shapefunction.dNdxi[0]=(eta - 1.0) * (1.0 - zeta) / 8.0;
            shapefunction.dNdxi[1]=(eta - 1.0) * (1.0 + zeta) / 8.0;
            shapefunction.dNdxi[2]=(1.0 - eta) * (1.0 + zeta) / 8.0;
            shapefunction.dNdxi[3]=(1.0 - eta) * (1.0 - zeta) / 8.0;
            shapefunction.dNdxi[4]=(1.0 + eta) * (1.0 - zeta) / 8.0;
            shapefunction.dNdxi[5]=(1.0 + eta) * (zeta - 1.0) / 8.0;
            shapefunction.dNdxi[6]=(-1.0 - eta) * (1.0 + zeta) / 8.0;
            shapefunction.dNdxi[7]=(1.0 + eta) * (1.0 + zeta) / 8.0;

            shapefunction.dNdeta[0]=(xi - 1.0) * (1.0 - zeta) / 8.0;
            shapefunction.dNdeta[1]=(xi - 1.0) * (1.0 + zeta) / 8.0;
            shapefunction.dNdeta[2]=(-1.0 - xi) * (1.0 + zeta) / 8.0;
            shapefunction.dNdeta[3]=(1.0 + xi) * (zeta - 1.0) / 8.0;
            shapefunction.dNdeta[4]=(1.0 + xi) * (1.0 - zeta) / 8.0;
            shapefunction.dNdeta[5]=(1.0 - xi) * (1.0 - zeta) / 8.0;
            shapefunction.dNdeta[6]=(1.0 - xi) * (1.0 + zeta) / 8.0;
            shapefunction.dNdeta[7]=(1.0 + xi) * (1.0 + zeta) / 8.0;

            shapefunction.dNdzeta[0]=(xi - 1.0) * (1.0 - eta) / 8.0;
            shapefunction.dNdzeta[1]=(1.0 - xi) * (1.0 - eta) / 8.0;
            shapefunction.dNdzeta[2]=(1.0 + xi) * (1.0 - eta) / 8.0;
            shapefunction.dNdzeta[3]=(1.0 + xi) * (eta - 1.0) / 8.0;
            shapefunction.dNdzeta[4]=(-1.0 - xi) * (1.0 + eta) / 8.0;
            shapefunction.dNdzeta[5]=(xi - 1.0) * (1.0 + eta) / 8.0;
            shapefunction.dNdzeta[6]=(1.0 - xi) * (1.0 + eta) / 8.0;
            shapefunction.dNdzeta[7]=(1.0 + xi) * (1.0 + eta) / 8.0;

            return shapefunction ;
        }

        struct JacobianData {
            std::array<std::array<double, 3>, 3> JacobianMatrix{};

            double detJ=0.0;

            std::array<std::array<double, 3>, 8> dNdX{};
        };

        [[nodiscard]]
        JacobianData computeH8Jacobian(
            const H8NodalValues& referenceCoordinates,
            const H8ShapeFunctions& shape
        ) {
            JacobianData J{};

            for(std::size_t i = 0;i < kH8NodeCount ; ++i) {
                J.JacobianMatrix[0][0] += shape.dNdxi[i] * referenceCoordinates[i][0];
                J.JacobianMatrix[0][1] += shape.dNdeta[i] * referenceCoordinates[i][0];
                J.JacobianMatrix[0][2] += shape.dNdzeta[i] * referenceCoordinates[i][0];
                J.JacobianMatrix[1][0] += shape.dNdxi[i] * referenceCoordinates[i][1];
                J.JacobianMatrix[1][1] += shape.dNdeta[i] * referenceCoordinates[i][1];
                J.JacobianMatrix[1][2] += shape.dNdzeta[i] * referenceCoordinates[i][1];
                J.JacobianMatrix[2][0] += shape.dNdxi[i] * referenceCoordinates[i][2];
                J.JacobianMatrix[2][1] += shape.dNdeta[i] * referenceCoordinates[i][2];
                J.JacobianMatrix[2][2] += shape.dNdzeta[i] * referenceCoordinates[i][2];
            }

            J.detJ = 
                J.JacobianMatrix[0][0] * J.JacobianMatrix[1][1] * J.JacobianMatrix[2][2] + J.JacobianMatrix[0][1] * J.JacobianMatrix[1][2] * J.JacobianMatrix[2][0] + J.JacobianMatrix[1][0] * J.JacobianMatrix[2][1] * J.JacobianMatrix[0][2] -
                J.JacobianMatrix[0][2] * J.JacobianMatrix[1][1] * J.JacobianMatrix[2][0] - J.JacobianMatrix[0][1] * J.JacobianMatrix[1][0] * J.JacobianMatrix[2][2] - J.JacobianMatrix[0][0] * J.JacobianMatrix[1][2] * J.JacobianMatrix[2][1];

            double scale = 0.0;
            for(std::size_t i =0 ;i < 3;++i) for(std::size_t j = 0;j < 3;++j) scale = std::max(scale,std::abs(J.JacobianMatrix[i][j]));
            const double tolerance = 100.0 * std::numeric_limits<double>::epsilon() * scale * scale *scale;

            if (!std::isfinite(J.detJ) || J.detJ <= tolerance) {
                throw std::invalid_argument(
                    "Invalid H8 Jacobian."
                );
            }

            std::array<std::array<double , 3>, 3> invJ{};
            invJ[0][0] = (J.JacobianMatrix[1][1] * J.JacobianMatrix[2][2] - J.JacobianMatrix[2][1] * J.JacobianMatrix[1][2])/J.detJ;
            invJ[1][0] = (J.JacobianMatrix[2][0] * J.JacobianMatrix[1][2] - J.JacobianMatrix[1][0] * J.JacobianMatrix[2][2])/J.detJ;
            invJ[2][0] = (J.JacobianMatrix[1][0] * J.JacobianMatrix[2][1] - J.JacobianMatrix[2][0] * J.JacobianMatrix[1][1])/J.detJ;
            invJ[0][1] = (J.JacobianMatrix[2][1] * J.JacobianMatrix[0][2] - J.JacobianMatrix[0][1] * J.JacobianMatrix[2][2])/J.detJ;
            invJ[1][1] = (J.JacobianMatrix[0][0] * J.JacobianMatrix[2][2] - J.JacobianMatrix[2][0] * J.JacobianMatrix[0][2])/J.detJ;
            invJ[2][1] = (J.JacobianMatrix[2][0] * J.JacobianMatrix[0][1] - J.JacobianMatrix[0][0] * J.JacobianMatrix[2][1])/J.detJ;
            invJ[0][2] = (J.JacobianMatrix[0][1] * J.JacobianMatrix[1][2] - J.JacobianMatrix[1][1] * J.JacobianMatrix[0][2])/J.detJ;
            invJ[1][2] = (J.JacobianMatrix[1][0] * J.JacobianMatrix[0][2] - J.JacobianMatrix[0][0] * J.JacobianMatrix[1][2])/J.detJ;
            invJ[2][2] = (J.JacobianMatrix[0][0] * J.JacobianMatrix[1][1] - J.JacobianMatrix[1][0] * J.JacobianMatrix[0][1])/J.detJ;


            for(std::size_t i = 0;i < kH8NodeCount ; ++i) {
                J.dNdX[i][0] = shape.dNdxi[i] * invJ[0][0] + shape.dNdeta[i] * invJ[1][0] + shape.dNdzeta[i] * invJ[2][0];
                J.dNdX[i][1] = shape.dNdxi[i] * invJ[0][1] + shape.dNdeta[i] * invJ[1][1] + shape.dNdzeta[i] * invJ[2][1];
                J.dNdX[i][2] = shape.dNdxi[i] * invJ[0][2] + shape.dNdeta[i] * invJ[1][2] + shape.dNdzeta[i] * invJ[2][2];
            }

            return J;

        }

        Matrix3D computeH8DisplacementGradient (
            const H8NodalValues& nodaldisplacements,
            const JacobianData& jacobian 
        )  {
            Matrix3D displacementgradient{};

            for(std::size_t i = 0 ;i < kH8NodeCount ; ++i) {
                displacementgradient[0][0] += jacobian.dNdX[i][0] * nodaldisplacements[i][0]; 
                displacementgradient[0][1] += jacobian.dNdX[i][1] * nodaldisplacements[i][0]; 
                displacementgradient[0][2] += jacobian.dNdX[i][2] * nodaldisplacements[i][0];
                displacementgradient[1][0] += jacobian.dNdX[i][0] * nodaldisplacements[i][1]; 
                displacementgradient[1][1] += jacobian.dNdX[i][1] * nodaldisplacements[i][1]; 
                displacementgradient[1][2] += jacobian.dNdX[i][2] * nodaldisplacements[i][1];
                displacementgradient[2][0] += jacobian.dNdX[i][0] * nodaldisplacements[i][2]; 
                displacementgradient[2][1] += jacobian.dNdX[i][1] * nodaldisplacements[i][2]; 
                displacementgradient[2][2] += jacobian.dNdX[i][2] * nodaldisplacements[i][2];
            }

            return displacementgradient;
        }

        struct H8GaussPointState {
            double detJ = 0.0;
            Matrix3D deformationGradient{};
            GreenLagrangeStrain_3D strain{};
            SecondPKStress_3D stress{};
            std::array<std::array<double, 3>, 8> dNdX{};
        };

        [[nodiscard]]
        H8GaussPointState computeH8GaussPointState (
            double xi,
            double eta,
            double zeta,
            const H8NodalValues& referenceCoordinates,
            const H8NodalValues& nodaldisplacements,
            const MaterialParameters& material
        ) {
            const H8ShapeFunctions shape = computeH8Shapefunctions(xi,eta,zeta);

            const JacobianData jacobian = computeH8Jacobian(referenceCoordinates,shape);

            const Matrix3D displacementgradient = computeH8DisplacementGradient(nodaldisplacements , jacobian);

            const Matrix3D deformationgradient = computeDeformationGradient_3D(displacementgradient);

            const GreenLagrangeStrain_3D strain = computeGreenLagrangeStrain_3D(deformationgradient);

            const SecondPKStress_3D stress = computeSecondPKStress_3D(strain , material);

            H8GaussPointState state{};
            state.detJ = jacobian.detJ;
            state.deformationGradient = deformationgradient;
            state.strain = strain;
            state.stress = stress;
            state.dNdX = jacobian.dNdX;
            
            return state;
        }

        [[nodiscard]]
        Matrix3D computeFirstPKStress_3D (
            const Matrix3D& deformationgradient,
            const SecondPKStress_3D& stress
        ) {
            Matrix3D firstPKStress{};

            firstPKStress[0][0] = deformationgradient[0][0] * stress.sxx + deformationgradient[0][1] * stress.sxy + deformationgradient[0][2] * stress.sxz;
            firstPKStress[0][1] = deformationgradient[0][0] * stress.sxy + deformationgradient[0][1] * stress.syy + deformationgradient[0][2] * stress.syz;
            firstPKStress[0][2] = deformationgradient[0][0] * stress.sxz + deformationgradient[0][1] * stress.syz + deformationgradient[0][2] * stress.szz;
            firstPKStress[1][0] = deformationgradient[1][0] * stress.sxx + deformationgradient[1][1] * stress.sxy + deformationgradient[1][2] * stress.sxz;
            firstPKStress[1][1] = deformationgradient[1][0] * stress.sxy + deformationgradient[1][1] * stress.syy + deformationgradient[1][2] * stress.syz;
            firstPKStress[1][2] = deformationgradient[1][0] * stress.sxz + deformationgradient[1][1] * stress.syz + deformationgradient[1][2] * stress.szz;
            firstPKStress[2][0] = deformationgradient[2][0] * stress.sxx + deformationgradient[2][1] * stress.sxy + deformationgradient[2][2] * stress.sxz;
            firstPKStress[2][1] = deformationgradient[2][0] * stress.sxy + deformationgradient[2][1] * stress.syy + deformationgradient[2][2] * stress.syz;
            firstPKStress[2][2] = deformationgradient[2][0] * stress.sxz + deformationgradient[2][1] * stress.syz + deformationgradient[2][2] * stress.szz;

            return firstPKStress;
        }
    }




    H8InternalForce computeH8InternalForce(
        const H8NodalValues& referenceCoordinates,
        const H8NodalValues& nodalDisplacements,
        const MaterialParameters& material
    ){
        H8InternalForce internalForce{};

        constexpr double gauss = 0.57735026918962576451;

        const std::array<double , 2> gausspoints = {-gauss,gauss};

        for(std::size_t i = 0 ;i < 2; ++i) {
            for(std::size_t j = 0; j < 2; ++j) {
                for(std::size_t l = 0;l < 2;++l) {
                    const double xi = gausspoints[i] ,eta = gausspoints[j],zeta = gausspoints[l];

                    const H8GaussPointState state = computeH8GaussPointState(xi , eta ,zeta , referenceCoordinates , nodalDisplacements ,material);

                    const Matrix3D firstPKStress = computeFirstPKStress_3D(state.deformationGradient,state.stress);

                    const double integrationWeight = state.detJ;

                    for(std::size_t k = 0; k < kH8NodeCount ; ++k) {
                        double forcex = 0.0, forcey = 0.0,forcez = 0.0;

                        for(std::size_t h = 0;h < 3 ; ++h) {
                            forcex += state.dNdX[k][h] * firstPKStress[0][h];
                            forcey += state.dNdX[k][h] * firstPKStress[1][h];
                            forcez += state.dNdX[k][h] * firstPKStress[2][h];
                        } 

                        internalForce[3*k] += forcex * integrationWeight;
                        internalForce[3*k+1] += forcey * integrationWeight;
                        internalForce[3*k+2] += forcez * integrationWeight;
                    }
                }
            }
        }

        return internalForce;
    }

    H8TangentStiffness computeH8TangentStiffness(
        const H8NodalValues& referenceCoordinates,
        const H8NodalValues& nodalDisplacements,
        const MaterialParameters& material
    ){
        H8TangentStiffness K{} ;

        const MaterialTangentMatrix_3D C = computeMaterialTangent_3D(material);

        constexpr double gauss = 0.57735026918962576451;

        const std::array<double , 2> gausspoints = {-gauss,gauss};

        std::array<std::array<double,kH8DofCount>,6> B{};

        for(std::size_t i = 0; i < 2; ++i) {
            for(std::size_t j = 0 ;j < 2; ++j) {
                for(std::size_t l = 0;l < 2;++l) {
                    const double xi = gausspoints[i],eta = gausspoints[j],zeta = gausspoints[l];

                    H8GaussPointState state = computeH8GaussPointState(xi,eta,zeta,referenceCoordinates,nodalDisplacements,material);

                    for(std::size_t k = 0 ;k < kH8NodeCount ;++k) {
                        const double Nx = state.dNdX[k][0],Ny = state.dNdX[k][1],Nz = state.dNdX[k][2]; 
                        for(std::size_t h = 0;h < 3;++h) {
                            size_t pos = 3 * k + h;
                            B[0][pos] = state.deformationGradient[h][0] * Nx;
                            B[1][pos] = state.deformationGradient[h][1] * Ny;
                            B[2][pos] = state.deformationGradient[h][2] * Nz;

                            B[3][pos] = state.deformationGradient[h][0] * Ny + state.deformationGradient[h][1] * Nx;
                            B[4][pos] = state.deformationGradient[h][0] * Nz + state.deformationGradient[h][2] * Nx;
                            B[5][pos] = state.deformationGradient[h][1] * Nz + state.deformationGradient[h][2] * Ny;
                        }
                    }

                    const double integrationWeight = state.detJ;


                    for(std::size_t k = 0 ;k < kH8DofCount ;++k) {
                        for(std::size_t h = 0; h < kH8DofCount ;++h) {
                            for(size_t r = 0;r < 6 ;++r) {
                                for(size_t s = 0 ;s < 6;++s) {
                                    K[k][h] += B[r][k] * C[r][s] * B[s][h] *integrationWeight;
                                }
                            }
                        }
                    }

                    const std::array<std::array<double,3> , 3> S={{
                        {state.stress.sxx , state.stress.sxy , state.stress.sxz} ,
                        {state.stress.sxy , state.stress.syy , state.stress.syz} ,
                        {state.stress.sxz , state.stress.syz , state.stress.szz}
                    }};

                    for(std::size_t k = 0 ;k < kH8NodeCount ;++k) {
                        for(std::size_t h = 0; h < kH8NodeCount ;++h) {
                            for(size_t r = 0;r < 3 ;++r) {
                                for(size_t s = 0 ;s < 3;++s) {
                                    K[k * 3][h * 3] += state.dNdX[k][r] * S[r][s] * state.dNdX[h][s] * integrationWeight;
                                    K[k * 3 +1][h * 3 + 1] += state.dNdX[k][r] * S[r][s] * state.dNdX[h][s] * integrationWeight;
                                    K[k * 3 + 2][h * 3 + 2] += state.dNdX[k][r] * S[r][s] * state.dNdX[h][s] * integrationWeight;
                                }
                            }
                        }
                    }
                }
            }
        }
        return K;
    }
    
    H8LumpedMass computeH8LumpedMass (
        const H8NodalValues& referenceCoordinates,
        double referenceDensity
    ) {
        if (!std::isfinite(referenceDensity) || referenceDensity <= 0.0) {
            throw std::invalid_argument(
                "Reference density must be finite and positive."
            );
        }


        H8LumpedMass mass{};

        constexpr double gauss = 0.57735026918962576451;

        const std::array<double , 2> gausspoints = {-gauss,gauss};

        for(std::size_t i = 0;i < 2;++i) {
            for(std::size_t j=0;j < 2;++j) {
                for(std::size_t l = 0 ;l < 2;++l) {
                    const double xi =gausspoints[i] ,eta = gausspoints[j],zeta = gausspoints[l];

                    const H8ShapeFunctions shape = computeH8Shapefunctions(xi,eta,zeta);

                    const JacobianData jacobian = computeH8Jacobian(referenceCoordinates,shape);

                    double m = referenceDensity * jacobian.detJ;

                    for(std::size_t k = 0;k < kH8NodeCount;++k) {
                        mass[3 * k] += m * shape.N[k];
                        mass[3 * k + 1] += m * shape.N[k];
                        mass[3 * k + 2] += m * shape.N[k];
                    }

                }
            }
        } 
        return mass;
    }

    double computeH8StrainEnergy (
        const H8NodalValues& referenceCoordinates,
        const H8NodalValues& nodalDisplacements,
        const MaterialParameters& material
    ) {
        double energy = 0.0;

        constexpr double gauss = 0.57735026918962576451;

        const std::array<double , 2> gausspoints = {-gauss,gauss};

        for(std::size_t i = 0;i < 2;++i) {
            for(std::size_t j = 0;j < 2;++j) {
                for(std::size_t l = 0;l < 2;++l) {
                    const double xi =gausspoints[i] ,eta = gausspoints[j],zeta = gausspoints[l];

                    H8GaussPointState state = computeH8GaussPointState(xi,eta,zeta,referenceCoordinates,nodalDisplacements,material);

                    energy += state.detJ * (state.stress.sxx * state.strain.exx + state.stress.syy * state.strain.eyy +state.stress.szz *state.strain.ezz + state.stress.sxy * state.strain.exy + state.stress.sxz * state.strain.exz + state.stress.syz * state.strain.eyz) / 2;
                }

            }
        } 
        if (!std::isfinite(energy) || energy < 0.0) {
            throw std::invalid_argument(
                "Energy must be finite and positive."
            );
        }

        return energy;
    }
}