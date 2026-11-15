#include "nonlinearfem/DynamicAnalysis.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>
#include <string>
#include <limits>

namespace NonlinearFEM {
    namespace {
        DynamicRecord SaveRecord(
            const ExplicitState& state,
            const GlobalLumpedMass& mass,
            const Q4Mesh& mesh,
            const MaterialParameters& material,
            double thickness
        ) {
            DynamicRecord record{};

            record.state = state;

            const size_t N = state.velocity.size();
            
            for(std::size_t i = 0;i < N;++i) {
                record.kineticEnergy += state.velocity[i] * state.velocity[i] * mass[i] / 2;
            }

            record.strainEnergy = assembleGlobalEnergy(mesh,state.displacement,material,thickness);

            if (!std::isfinite(record.kineticEnergy) || record.kineticEnergy < 0.0 || !std::isfinite(record.strainEnergy) || record.strainEnergy < 0.0) {
                throw std::runtime_error(
                    "Dynamic energies must be finite and nonnegative."
                );
            }

            return record;
        }
    }



    [[nodiscard]] 
    DynamicAnalysisResult solveDynamicAnalysis(
        const Q4Mesh& mesh,
        const MaterialParameters& material,
        double referenceDensity,
        double thickness ,
        const std::vector<bool>& constrained, 
        const ExplicitState& initialState,
        const ExternalForceFunction& externalForce,
        const DynamicAnalysisOptions& options
    ) {
        if (!std::isfinite(options.timeStep) || options.timeStep <= 0.0) {
            throw std::invalid_argument(
                "Time step must be finite and positive."
            );
        }

        if (options.numberOfSteps == 0 || options.outputEvery == 0) {
            throw std::invalid_argument(
                "Number of steps and output interval must be positive."
            );
        }

        if (!std::isfinite(referenceDensity) ||  referenceDensity <= 0.0) {
            throw std::invalid_argument(
                "Reference density must be finite and positive."
            );
        }

        if (!std::isfinite(thickness) || thickness <= 0.0) {
            throw std::invalid_argument(
                "Thickness must be finite and positive."
            );
        }

        if (!std::isfinite(initialState.time)) {
            throw std::invalid_argument(
                "Initial time must be finite."
            );
        }

        if (!externalForce) {
            throw std::invalid_argument(
                "External force function must be valid."
            );
        }

        if (mesh.nodes.size() > std::numeric_limits<std::size_t>::max() / 2) {
            throw std::overflow_error(
                "Global DOF count overflow."
            );
        }

        const std::size_t N = 2 * mesh.nodes.size();

        if (N == 0 || initialState.displacement.size() != N || initialState.velocity.size() != N || constrained.size() != N) {
            throw std::invalid_argument(
                "Initial state or constraints do not match mesh DOFs."
            );
        }

        const double duration = options.timeStep * static_cast<double>(options.numberOfSteps);

        const double endTime = initialState.time + duration;

        if (!std::isfinite(duration) || !std::isfinite(endTime) || endTime <= initialState.time) {
            throw std::invalid_argument(
                "Analysis duration must produce a finite, increasing time."
            );
        }

        const GlobalLumpedMass mass = assembleGlobalLumpedMass(mesh,referenceDensity,thickness);

        ExplicitState state = initialState;

        const InternalForceFunction internalforce= [&mesh,&material,thickness](const DynamicVector& u) {
            return assembleInternalForce(mesh,u,material,thickness);
        };

        InitializeExplicitAcceleration(state,mass,constrained,internalforce,externalForce);

        DynamicAnalysisResult result{};

        result.history.push_back(SaveRecord(state,mass,mesh,material,thickness));

        for(std::size_t step = 0;step < options.numberOfSteps ;++step) {
            advanceCentralDifference(state,mass,constrained,options.timeStep,internalforce,externalForce);

            if((step + 1) % options.outputEvery == 0 || (step + 1) == options.numberOfSteps) result.history.push_back(SaveRecord(state,mass,mesh,material,thickness));
        }

        result.finalState = std::move(state);

        return result;
    }
}