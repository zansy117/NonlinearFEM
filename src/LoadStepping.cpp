#include "nonlinearfem/LoadStepping.h"

#include <cstddef>
#include <string>
#include <vector>
#include <algorithm>

namespace NonlinearFEM {
    [[nodiscard]]
    LoadSteppingResult solveLoadStepping (
        const Q4Mesh& mesh,
        const MaterialParameters& material,
        const NodalLoads& loads,
        const DisplacementConstraints& constraints,
        const GlobalDisplacement& initialdisplacement,
        const LoadSteppingOptions& options,
        double thickness 
    ) {
        LoadSteppingResult result{};

        result.displacement = initialdisplacement;

        if (options.initialIncrement <= 0.0 ||
            options.minimumIncrement <= 0.0 ||
            options.maximumIncrement <= 0.0 ||
            options.cutbackFactor <= 0.0 ||
            options.cutbackFactor >= 1.0 ||
            options.growthFactor < 1.0 ||
            options.minimumIncrement > options.maximumIncrement ||
            options.initialIncrement < options.minimumIncrement ||
            options.initialIncrement > options.maximumIncrement ||
            options.maxLoadSteps == 0) {

            result.failureReason =
                "Invalid load stepping options.";

            return result;
        }

        double loadFactor = 0.0;

        double increment = options.initialIncrement;

        GlobalDisplacement displacementNow = initialdisplacement;

        for(std::size_t step = 0;step < options.maxLoadSteps ;++step) {
            std::size_t retries = 0;

            bool stepConverged = false;

            LoadStepRecord record{};

            while(retries <= options.maxRetries) {
                const double targetLoadFactor = std::min(1.0, loadFactor + increment);
                const double actualIncrement = targetLoadFactor - loadFactor;

                record.loadFactor = targetLoadFactor;
                record.loadIncrement = actualIncrement;

                NodalLoads scaleloads{};
                scaleloads.resize(loads.size());

                for(size_t i = 0;i < loads.size();++i) {
                    scaleloads[i] = loads[i];
                    scaleloads[i].value *= targetLoadFactor;
                }

                DisplacementConstraints scaleconstraints{};
                scaleconstraints.resize(constraints.size());

                for(size_t i = 0;i < constraints.size();++i) {
                    scaleconstraints[i] = constraints[i];
                    scaleconstraints[i].value *= targetLoadFactor;
                }

                const GlobalInternalForce externalForce = assembleExternalForce(mesh.nodes.size(),scaleloads);

                NewtonRaphsonResult newtonResult = solveNewtonRaphson(mesh,material,externalForce,scaleconstraints,displacementNow,options.newtonOptions,thickness);

                if(newtonResult.converged) {
                    stepConverged = true;

                    displacementNow = newtonResult.displacement;

                    loadFactor = targetLoadFactor;

                    record.converged = true;
                    record.iterations = newtonResult.iterations;
                    record.residualNorm = newtonResult.residualNorm;
                    record.retries = retries;

                    break;
                }    

                record.failureReason = newtonResult.failureReason;
                
                retries++;
                if (retries > options.maxRetries) break;

                increment *= options.cutbackFactor;
                if (increment < options.minimumIncrement) break;
            }    

                if (!stepConverged) {
                    record.converged = false;

                    record.retries = retries;

                    result.history.push_back(record);

                    result.displacement = displacementNow;

                    result.finalLoadFactor = loadFactor;

                    result.failureReason = "Load increment failed: " + record.failureReason;

                    return result;
                }

                result.history.push_back(record);
                result.displacement = displacementNow;
                result.finalLoadFactor = loadFactor;

                if (loadFactor >= 1.0 - 1.0e-12) {
                    result.converged = true;
                    return result;
                }

                increment = std::min(options.maximumIncrement,increment * options.growthFactor);
        }
        result.displacement = displacementNow;
        result.finalLoadFactor = loadFactor;
        result.failureReason = "Maximum number of load steps reached.";
        return result;
    }
    
}