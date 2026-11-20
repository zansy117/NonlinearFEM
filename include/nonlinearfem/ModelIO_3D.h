#pragma once
#include "nonlinearfem/LoadStepping_3D.h"
#include "nonlinearfem/DynamicAnalysis_3D.h"
#include <filesystem>
#include <istream>
namespace NonlinearFEM {
struct AnalysisModel_3D {
    H8Mesh mesh;
    MaterialParameters material;
    NodalLoads_3D loads;
    DisplacementConstraints_3D constraints;
    LoadSteppingOptions options;
};
AnalysisModel_3D readModel_3D(std::istream& input);
AnalysisModel_3D readModelFile_3D(const std::filesystem::path& path);
void writeDynamicResults_3D(const H8Mesh& mesh,const DynamicAnalysisResult& result,
    std::size_t observationNode,const std::filesystem::path& directory);
void writeResults_3D(const AnalysisModel_3D& model,const LoadSteppingResult_3D& result,
    const std::filesystem::path& directory);
void runDynamicModel_3D(const AnalysisModel_3D& model,const std::filesystem::path& configuration,
    const std::filesystem::path& output);
}
