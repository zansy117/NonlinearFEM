#pragma once
#include "nonlinearfem/LoadStepping.h"
#include <filesystem>
#include <istream>

namespace NonlinearFEM {
struct AnalysisModel {
    Q4Mesh mesh;
    MaterialParameters material;
    double thickness = 1.0;
    NodalLoads loads;
    DisplacementConstraints constraints;
    LoadSteppingOptions options;
};
// Strict, line-oriented NFEM 1 format; node indices are zero-based.
AnalysisModel readModel(std::istream& input);
AnalysisModel readModelFile(const std::filesystem::path& path);
// Only fully converged results are exported. Refuses to overwrite an existing run.
void writeResults(const AnalysisModel& model, const LoadSteppingResult& result,
                  const std::filesystem::path& directory);
}
