#pragma once

#include <functional>
#include <vector>

namespace NonlinearFEM {
    using DynamicVector = std::vector<double>;

    struct ExplicitState {
        double time = 0.0;

        DynamicVector displacement;
        DynamicVector velocity;
        DynamicVector acceleration;  //加速度
    };

    using InternalForceFunction = std::function< DynamicVector(const DynamicVector&) >;
    using ExternalForceFunction = std::function< DynamicVector(double) >;

    void InitializeExplicitAcceleration (
        ExplicitState& state,
        const DynamicVector& mass,
        const std::vector<bool>& constrained,
        const InternalForceFunction& internalForce,
        const ExternalForceFunction& externalForce
    );

    void advanceCentralDifference (
        ExplicitState& state,
        const DynamicVector& mass,
        const std::vector<bool>& constrained,
        double dt,
        const InternalForceFunction& internalForce,
        const ExternalForceFunction& externalForce
    );
}