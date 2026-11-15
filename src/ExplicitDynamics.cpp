#include "nonlinearfem/ExplicitDynamics.h"

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
} 