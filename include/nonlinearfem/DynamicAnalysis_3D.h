#pragma once

#include "nonlinearfem/Assembly_3D.h"
#include "nonlinearfem/DynamicAnalysis.h"

namespace NonlinearFEM {

    [[nodiscard]] 
    DynamicAnalysisResult solveDynamicAnalysis_3D(
        const H8Mesh& mesh,
        const MaterialParameters& material,
        double referenceDensity,
        const std::vector<bool>& constrained, 
        const ExplicitState& initialState,
        const ExternalForceFunction& externalForce,
        const DynamicAnalysisOptions& options
    );
}