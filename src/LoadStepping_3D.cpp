#include "nonlinearfem/LoadStepping_3D.h"

#include <cstddef>
#include <string>
#include <vector>
#include <algorithm>
#include <cmath>

namespace NonlinearFEM {
    [[nodiscard]]
    LoadSteppingResult_3D solveLoadStepping_3D (
        const H8Mesh& mesh,
        const MaterialParameters& material,
        const NodalLoads_3D& loads,
        const DisplacementConstraints_3D& constraints,
        const GlobalDisplacement_3D& initialdisplacement,
        const LoadSteppingOptions& options
    ) {
        LoadSteppingResult_3D result{};

        result.displacement = initialdisplacement;

        if (!std::isfinite(options.initialIncrement) ||
            !std::isfinite(options.minimumIncrement) ||
            !std::isfinite(options.maximumIncrement) ||
            !std::isfinite(options.cutbackFactor) ||
            !std::isfinite(options.growthFactor) ||
            options.initialIncrement <= 0.0 ||
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

        GlobalDisplacement_3D displacementNow = initialdisplacement;

        for(std::size_t step = 0;step < options.maxLoadSteps ;++step) {
            std::size_t retries = 0;

            bool stepConverged = false;

            LoadStepRecord record{};

            while(retries <= options.maxRetries) {
                const double targetLoadFactor = std::min(1.0, loadFactor + increment);
                const double actualIncrement = targetLoadFactor - loadFactor;

                record.loadFactor = targetLoadFactor;
                record.loadIncrement = actualIncrement;

                if (!std::isfinite(targetLoadFactor) ||
                    !std::isfinite(actualIncrement) || actualIncrement <= 0.0) {
                    record.failureReason = "Load factor cannot advance.";
                    break;
                }

                NodalLoads_3D scaleloads{};
                scaleloads.resize(loads.size());

                for(size_t i = 0;i < loads.size();++i) {
                    scaleloads[i] = loads[i];
                    scaleloads[i].value *= targetLoadFactor;
                }

                DisplacementConstraints_3D scaleconstraints{};
                scaleconstraints.resize(constraints.size());

                for(size_t i = 0;i < constraints.size();++i) {
                    scaleconstraints[i] = constraints[i];
                    scaleconstraints[i].value *= targetLoadFactor;
                }

                const GlobalInternalForce_3D externalForce = assembleExternalForce_3D(mesh.nodes.size(),scaleloads);

                NewtonRaphsonResult_3D newtonResult = solveNewtonRaphson_3D(mesh,material,externalForce,scaleconstraints,displacementNow,options.newtonOptions);

                record.iterations = newtonResult.iterations;
                record.residualNorm = newtonResult.residualNorm;

                if(newtonResult.converged) {
                    stepConverged = true;

                    displacementNow = newtonResult.displacement;

                    loadFactor = targetLoadFactor;

                    record.converged = true;
                    record.failureReason.clear();
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