#pragma once
#include "nonlinearfem/DynamicAnalysis.h"
#include <filesystem>
namespace NonlinearFEM {
// Creates a new directory; refuses to overwrite an existing run.
void writeDynamicResults(const Q4Mesh& mesh,
    const DynamicAnalysisResult& result, std::size_t observationNode,
    const std::filesystem::path& directory);
}
