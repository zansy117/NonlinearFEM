#include"nonlinearfem/Solver_3D.h"

#include<cmath>
#include<vector>
#include<string>
#include <algorithm>
#include <cstddef>
#include <utility>
#include <limits>



namespace NonlinearFEM {
    namespace {
        [[nodiscard]]
        bool GaussElimination (
            GlobalTangentStiffness_3D K,
            GlobalInternalForce_3D F,
            GlobalDisplacement_3D& u 
        ) {
            std::size_t N = u.size();

            if (N == 0 || K.size() != N || F.size() != N) {
                return false;
            }

            for(std::size_t i = 0;i < N;++i) {
                if (K[i].size() != N || !std::isfinite(F[i])) {
                    return false;
                }

                for(std::size_t j = 0;j < N;++j) {
                    if (!std::isfinite(K[i][j])) return false;
                }
            }

            for(std::size_t i = 0;i < N;++i) {
                std::size_t tem = i;

                for(std::size_t j = i + 1;j < N;++j) if(fabs(K[tem][i]) < fabs(K[j][i])) tem = j;

                if (!std::isfinite(K[tem][i]) || std::abs(K[tem][i]) < 1e-14) {
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

                    for(std::size_t k = i;k < N;++k) K[j][k] -= K[i][k] * div;

                    F[j] -= F[i]*div;
                }
            }

            u[N - 1] = F[N - 1];
            if (!std::isfinite(u[N - 1])) {
                return false;
            }

            for(std::size_t i = N - 1;i > 0;) {
                --i;
                u[i] = F[i];

                for(std::size_t j = i + 1;j < N;++j) u[i] -= K[i][j] * u[j];

                if (!std::isfinite(u[i])) {
                    return false;
                }
            }
            return true;
        }
    }





    NewtonRaphsonResult_3D solveNewtonRaphson_3D(
        const H8Mesh& mesh,
        const MaterialParameters& material,
        const GlobalInternalForce_3D& externalForce,
        const DisplacementConstraints_3D& constraints,
        const GlobalDisplacement_3D& initialDisplacement,
        const NewtonRaphsonOptions& options
    ) {
        NewtonRaphsonResult_3D result{};
        result.displacement = initialDisplacement;

        const size_t N = externalForce.size();

        if (mesh.nodes.size() > std::numeric_limits<std::size_t>::max() / 3) {
            result.failureReason = "Global DOF count overflow.";
            return result;
        }

        if (N == 0 || N != 3 * mesh.nodes.size() || initialDisplacement.size() != N) {
            result.failureReason = "Invalid system dimensions.";
            return result;
        }

        if (!std::isfinite(options.residualTolerance) || options.residualTolerance <= 0.0 ||
            !std::isfinite(options.displacementTolerance) || options.displacementTolerance <= 0.0) {
            result.failureReason = "Invalid Newton-Raphson tolerances.";
            return result;
        }

        for(std::size_t i = 0;i < N;++i) {
            if (!std::isfinite(externalForce[i]) || !std::isfinite(initialDisplacement[i])) {
                result.failureReason = "External force and initial displacement must be finite.";
                return result;
            }
        }

        std::vector<bool> constrained;
        constrained.assign(N,false);

        for(std:: size_t i = 0;i < constraints.size() ;++i) {
            if (constraints[i].node >= mesh.nodes.size() || constraints[i].direction >= 3) {
                result.failureReason = "Invalid displacement constraint position.";
                return result;
            }

            if (!std::isfinite(constraints[i].value)) {
                result.failureReason = "Prescribed displacement must be finite.";
                return result;
            }

            const std::size_t position = 3 * constraints[i].node + constraints[i].direction;

            if (constrained[position]) {
                result.failureReason = "Duplicate displacement constraint.";
                return result;
            }

            constrained[position] = true;
        }

        double externalForceNormSquared = 0.0;

        for(std::size_t i = 0;i < N;++i) {
            if(constrained[i]) continue;
            externalForceNormSquared += externalForce[i] * externalForce[i];
        }
        double externalForceNorm = std::sqrt(externalForceNormSquared);

        if (!std::isfinite(externalForceNorm)) {
            result.failureReason = "External force norm overflow.";
            return result;
        }

        for(std::size_t iteration = 0;iteration <= options.maxIterations ;++iteration) {
            GlobalAssemblyResult_3D assemblyresult = assembleGlobalSystem_3D(mesh,result.displacement,material);

            GlobalInternalForce_3D residual{};
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

            if (!std::isfinite(residualNormSquared) || !std::isfinite(internalForceNormSquared)) {
                result.failureReason = "Non-finite force or residual norm.";
                return result;
            }

            double internalForceNorm = std::sqrt(internalForceNormSquared);

            double residualNorm = std::sqrt(residualNormSquared) / std::max(std::max(1.0 , externalForceNorm) , internalForceNorm);

            result.residualNorm = residualNorm;
            result.iterations = iteration;

            double constraintError = 0.0;

            for(std::size_t i = 0; i < constraints.size() ;++i) {
                std::size_t position = constraints[i].node * 3 + constraints[i].direction;

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

            GlobalTangentStiffness_3D stiffness = assemblyresult.tangentStiffness;

            applyDisplacementConstraints_3D(stiffness,residual,result.displacement,constraints);

            GlobalDisplacement_3D U{};
            U.assign(N,0.0);

            if (!GaussElimination(stiffness,residual,U)) {
                result.failureReason = "Linear system solution failed.";
                return result;
            }

            double displacementIncrementNorm = 0.0;

            for(std::size_t i = 0;i < N;++i) {
                displacementIncrementNorm = std::hypot(displacementIncrementNorm,U[i]);
            }

            if (!std::isfinite(displacementIncrementNorm)) {
                result.failureReason = "Displacement increment norm overflow.";
                return result;
            }

            result.displacementIncrementNorm = displacementIncrementNorm;

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