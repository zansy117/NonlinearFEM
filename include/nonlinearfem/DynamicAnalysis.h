#pragma once

#include "nonlinearfem/Assembly.h"
#include "nonlinearfem/ExplicitDynamics.h"

namespace NonlinearFEM {
    struct DynamicAnalysisOptions {
        double timeStep = 1.0e-4;
        std::size_t numberOfSteps = 1000;
        std::size_t outputEvery = 10;
    };

    struct DynamicRecord {
        ExplicitState state;

        double kineticEnergy = 0.0;
        double strainEnergy = 0.0;
    };

    struct DynamicAnalysisResult {
        ExplicitState finalState;
        std::vector<DynamicRecord> history;
    };

    [[nodiscard]] 
    DynamicAnalysisResult solveDynamicAnalysis(
        const Q4Mesh& mesh,
        const MaterialParameters& material,
        double referenceDensity,
        double thickness,
        const std::vector<bool>& constrained, 
        const ExplicitState& initialState,
        const ExternalForceFunction& externalForce,
        const DynamicAnalysisOptions& options
    );
}