#include "nonlinearfem/ExplicitDynamics.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <string>

namespace NonlinearFEM {
    namespace {
        void validateFiniteVector(
            const DynamicVector& values,
            const char* name
        ) {
            for (std::size_t i = 0; i < values.size(); ++i) {
                if (!std::isfinite(values[i])) {
                    throw std::invalid_argument(
                        std::string(name) + "[" + std::to_string(i) +
                        "] must be finite."
                    );
                }
            }       
        }

        void validateExplicitState(
            const ExplicitState& state,
            const DynamicVector& mass,
            const std::vector<bool>& constrained,
            const InternalForceFunction& internalForce,
            const ExternalForceFunction& externalForce,
            bool requireAcceleration
        ) {
            const std::size_t N = state.displacement.size();

            if (N == 0 || state.velocity.size() != N || mass.size() != N || constrained.size() != N) {
                throw std::invalid_argument(
                    "Inconsistent or empty dynamic system dimensions."
                );
            }

            if (!std::isfinite(state.time)) {
                throw std::invalid_argument("Time must be finite.");
            }

            if (!internalForce || !externalForce) {
                throw std::invalid_argument(
                    "Force functions must be valid."
                );
            }

            validateFiniteVector(state.displacement, "Displacement");
            validateFiniteVector(state.velocity, "Velocity");
            validateFiniteVector(mass, "Mass");

            if (requireAcceleration) {
                if (state.acceleration.size() != N) {
                    throw std::invalid_argument(
                        "Acceleration size does not match system DOFs."
                    );
                }

                validateFiniteVector(state.acceleration, "Acceleration");
            }

            for (std::size_t i = 0; i < N; ++i) {
                if (mass[i] < 0.0 ||
                    (!constrained[i] && mass[i] == 0.0)) {
                    throw std::invalid_argument(
                        "Mass must be nonnegative; free DOFs require positive mass."
                    );
                }

                if (constrained[i]) {
                    if (state.displacement[i] != 0.0 || state.velocity[i] != 0.0 || (requireAcceleration && state.acceleration[i] != 0.0)) {
                        throw std::invalid_argument(
                            "Fixed DOFs must have zero displacement, "
                            "velocity and acceleration."
                        );
                    }
                }
            }
        }

        void validateForceVector(
            const DynamicVector& force,
            std::size_t expectedSize,
            const char* name
        ) {
            if (force.size() != expectedSize) {
                throw std::invalid_argument(
                    std::string(name) + " size does not match system DOFs."
                );
            }

            validateFiniteVector(force, name);
        }

        void validateTimeStep(double time, double dt) {
            if (!std::isfinite(dt) || dt <= 0.0) {
                throw std::invalid_argument(
                    "Time step must be finite and positive."
                );
            }

            const double nextTime = time + dt;

            if (!std::isfinite(time) || !std::isfinite(nextTime) || nextTime <= time) {
                throw std::invalid_argument(
                    "Time step must produce a finite, increasing time."
                );
            }
        }

    }

    void InitializeExplicitAcceleration (
        ExplicitState& state,
        const DynamicVector& mass,
        const std::vector<bool>& constrained,
        const InternalForceFunction& internalForce,
        const ExternalForceFunction& externalForce
    ){
        validateExplicitState(state, mass, constrained,internalForce, externalForce,false);

        const std::size_t N = mass.size(); 

        DynamicVector acceleration{};
        acceleration.assign(N,0.0);

        DynamicVector intF = internalForce(state.displacement);
        DynamicVector extF = externalForce(state.time);

        validateForceVector(intF,N,"Internal force");
        validateForceVector(extF,N,"External force");

        for(std::size_t i = 0;i < N;++i) {
            if(constrained[i]) continue;

            acceleration[i] = (extF[i] - intF[i]) / mass[i];
        }

        validateFiniteVector(acceleration, "Initial acceleration");

        state.acceleration = std::move(acceleration);
    }

    void advanceCentralDifference (
        ExplicitState& state,
        const DynamicVector& mass,
        const std::vector<bool>& constrained,
        double dt,
        const InternalForceFunction& internalForce,
        const ExternalForceFunction& externalForce
    ) {
        validateExplicitState(state, mass, constrained,internalForce, externalForce,true);

        validateTimeStep(state.time, dt);

        const std::size_t N = mass.size();

        ExplicitState next = state;
        next.time += dt;

        DynamicVector halfVelocity{}; //半步速度
        halfVelocity.assign(N,0.0);

        for(std::size_t i = 0;i < N;++i) {
            if(constrained[i]) {
                next.displacement[i] = 0.0;
                next.velocity[i] = 0.0;
                next.acceleration[i] = 0.0;
                continue;
            }

            halfVelocity[i] = state.velocity[i] + dt * state.acceleration[i] / 2;
            next.displacement[i] += halfVelocity[i] * dt;
        }
        validateFiniteVector(next.displacement, "New displacement");

        DynamicVector intF = internalForce(next.displacement);
        DynamicVector extF = externalForce(next.time);

        validateForceVector(intF, N, "Internal force");
        validateForceVector(extF, N, "External force");

        for(std::size_t i = 0;i < N;++i) {
            if(constrained[i]) continue;

            next.acceleration[i] = (extF[i] - intF[i]) / mass[i];

            next.velocity[i] = halfVelocity[i] + next.acceleration[i] * dt / 2;
        }

        validateFiniteVector(next.displacement, "New displacement");
        validateFiniteVector(next.velocity, "New velocity");
        validateFiniteVector(next.acceleration, "New acceleration");

        state = std::move(next);
    }
    [[nodiscard]]
    double estimateStableTimeStep(
        const std::vector<std::vector<double>>& stiffness,
        const DynamicVector& mass,
        const std::vector<bool>& constrained,
        double safetyFactor
    ) {
        const std::size_t N = constrained.size();

        if (N == 0 || mass.size() != N || stiffness.size() != N) {
            throw std::invalid_argument("Inconsistent or empty time-step system dimensions.");
        }
        if (!std::isfinite(safetyFactor) || safetyFactor <= 0.0 || safetyFactor >= 1.0) {
            throw std::invalid_argument("Safety factor must be finite and in (0, 1).");
        }
        bool hasFreeDof = false;
        for (std::size_t i = 0; i < N; ++i) {
            if (stiffness[i].size() != N) {
                throw std::invalid_argument("Stiffness matrix must be square.");
            }
            if (!std::isfinite(mass[i]) || mass[i] < 0.0 || (!constrained[i] && mass[i] == 0.0)) {
                throw std::invalid_argument("Mass must be finite and nonnegative; free DOFs require positive mass.");
            }
            if (!constrained[i]) hasFreeDof = true;
            for (double value : stiffness[i]) {
                if (!std::isfinite(value)) {
                    throw std::invalid_argument("Stiffness entries must be finite.");
                }
            }
        }
        if (!hasFreeDof) {
            throw std::invalid_argument("At least one free DOF is required.");
        }
        // Compare normalized entries to avoid overflow in the symmetry check.
        for (std::size_t i = 0; i < N; ++i) {
            for (std::size_t j = i + 1; j < N; ++j) {
                const double scale = std::max(1.0, std::max(std::abs(stiffness[i][j]), std::abs(stiffness[j][i])));
                if (std::abs(stiffness[i][j] / scale - stiffness[j][i] / scale) > 1.0e-12) {
                    throw std::invalid_argument("Stiffness matrix must be symmetric.");
                }
            }
        }

        double answer = 0.0;
        for(size_t i = 0;i < N;++i) {
            if(constrained[i]) continue;

            double temp = 0.0;

            for(size_t j = 0;j < N;++j) {
                if(constrained[j]) continue;

                const double denominator = std::sqrt(mass[i]) * std::sqrt(mass[j]);
                if (!std::isfinite(denominator) || denominator <= 0.0) {
                    throw std::overflow_error("Invalid mass normalization denominator.");
                }
                double res = stiffness[i][j] / denominator;
                if (!std::isfinite(res)) {
                    throw std::overflow_error("Non-finite mass-normalized stiffness.");
                }

                temp += std::abs(res);
                if (!std::isfinite(temp)) {
                    throw std::overflow_error("Stiffness row sum overflow.");
                }
            }

            answer = std::max(answer,temp);
        }
        
        if (!std::isfinite(answer) || answer <= 0.0) {
            throw std::invalid_argument("Eigenvalue bound must be finite and positive.");
        }

        double dt = safetyFactor * 2 / std::sqrt(answer);
        if (!std::isfinite(dt) || dt <= 0.0) {
            throw std::overflow_error("Estimated time step must be finite and positive.");
        }

        return dt;
    }
} 