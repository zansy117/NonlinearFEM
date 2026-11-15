#pragma once
#include "nonlinearfem/ModelIO.h"
namespace NonlinearFEM {
void runDynamicModel(const AnalysisModel& model,
    const std::filesystem::path& configuration,
    const std::filesystem::path& output);
}
