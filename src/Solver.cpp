#include"nonlinearfem/Solver.h"

#include<cmath>
#include<vector>
#include<string>
#include <algorithm>
#include <cstddef>
#include <utility>



namespace NonlinearFEM {
    namespace {
        [[nodiscard]]
        bool GaussElimination (
            GlobalTangentStiffness K,
            GlobalInternalForce F,
            GlobalDisplacement& u 
        ) {
            std::size_t N = u.size();

            if (N == 0 || K.size() != N || F.size() != N) {
                return false;
            }

            for(std::size_t i = 0;i < N;++i) {
                std::size_t tem = i;

                for(std::size_t j = i + 1;j < N;++j) if(fabs(K[tem][i]) < fabs(K[j][i])) tem = j;

                if (std::abs(K[tem][i]) < 1e-14) {
                    return false;
                }

                if(i != tem) {
                    std::swap(K[i],K[tem]);
                    std::swap(F[i],F[tem]);
                }

                double div = K[i][i];


                for(size_t j = i;j < N;++j) {
                    K[i][j] /= div;
                }

                F[i] /= div;

                for(size_t j = i + 1;j < N;++j) {
                    double div = K[j][i];

                    for(int k = i;k < N;++k) K[j][k] -= K[i][k] * div;


                    F[j] -= F[i]*div;
                }
            }

            u[N - 1] = F[N - 1];
            if (!std::isfinite(u[N - 1])) {
                return false;
            }

            for(int i = N - 2;i >= 0;--i) {
                u[i] = F[i];

                for(int j = i + 1;j < N;++j) u[i] -= K[i][j] * u[j];

                if (!std::isfinite(u[i])) {
                    return false;
                }
            }
            return true;
        }
    }





    NewtonRaphsonResult solveNewtonRaphson(
        const Q4Mesh& mesh,
        const MaterialParameters& material,
        const GlobalInternalForce& externalForce,
        const DisplacementConstraints& constraints,
        const GlobalDisplacement& initialDisplacement,
        const NewtonRaphsonOptions& options,
        double thickness
    ) {
        NewtonRaphsonResult result{};
        result.displacement = initialDisplacement;

        const size_t N = externalForce.size();

        if (N == 0 || N != 2 * mesh.nodes.size() || initialDisplacement.size() != N) {
            result.failureReason = "Invalid system dimensions.";
            return result;
        }

        std::vector<bool> constrained;
        constrained.assign(N,false);

        for(std:: size_t i = 0;i < constraints.size() ;++i) {
            constrained[2 * constraints[i].node + constraints[i].direction] = true ;
        }

        double externalForceNormSquared = 0.0;

        for(std::size_t i = 0;i < N;++i) {
            if(constrained[i]) continue;
            externalForceNormSquared += externalForce[i] * externalForce[i];
        }
        double externalForceNorm = std::sqrt(externalForceNormSquared);

        for(std::size_t iteration = 0;iteration <= options.maxIterations ;++iteration) {
            GlobalAssemblyResult assemblyresult = assembleGlobalSystem(mesh,result.displacement,material,thickness);

            GlobalInternalForce residual{};
            residual.assign(N, 0.0);

            for(std::size_t i = 0 ;i < N ;++i) {
                residual[i] = externalForce[i] - assemblyresult.internalForce[i];
            }


            double residualNormSquared = 0.0;
            double internalForceNormSquared = 0.0;

            for(std::size_t i = 0;i < N;++i) {
                internalForceNormSquared += assemblyresult.internalForce[i] * assemblyresult.internalForce[i];
                if(constrained[i]) continue;
                residualNormSquared += residual[i] * residual[i];
            }

            double internalForceNorm = std::sqrt(internalForceNormSquared);

            double residualNorm = std::sqrt(residualNormSquared) / std::max(std::max(1.0 , externalForceNorm) , internalForceNorm);

            result.residualNorm = residualNorm;
            result.iterations = iteration;

            double constraintError = 0.0;

            for(std::size_t i = 0; i < constraints.size() ;++i) {
                std::size_t position = constraints[i].node * 2 + constraints[i].direction;

                constraintError = std::max(constraintError,std::abs(result.displacement[position] - constraints[i].value));
            }

            if (residualNorm <= options.residualTolerance && constraintError <= options.displacementTolerance) {
                result.converged = true;
                result.iterations = iteration;
                return result;
            }

            if (iteration == options.maxIterations) {
                result.failureReason = "Maximum Newton-Raphson iterations reached.";
                return result;
            }

            GlobalTangentStiffness stiffness = assemblyresult.tangentStiffness;

            applyDisplacementConstraints(stiffness,residual,result.displacement,constraints);

            GlobalDisplacement U{};
            U.assign(N,0.0);

            if (!GaussElimination(stiffness,residual,U)) {
                result.failureReason = "Linear system solution failed.";
                return result;
            }

            for (std::size_t i = 0; i < N; ++i) {
                result.displacement[i] += U[i];

                if (!std::isfinite(result.displacement[i])) {
                    result.failureReason = "Non-finite displacement detected.";
                    return result;
                }
            }
        }
        return result;
    }
}