#include "nonlinearfem/DynamicOutput.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdexcept>
namespace NonlinearFEM {
namespace {
std::ofstream open(const std::filesystem::path& path) {
    std::ofstream file(path);
    file.exceptions(std::ios::failbit | std::ios::badbit);
    file << std::setprecision(17);
    return file;
}
}
void writeDynamicResults(const Q4Mesh& mesh,
    const DynamicAnalysisResult& result, std::size_t observationNode,
    const std::filesystem::path& directory) {
    if (mesh.nodes.empty() || mesh.elements.empty() || result.history.empty() ||
        observationNode >= mesh.nodes.size())
        throw std::invalid_argument("Invalid mesh, observation node or empty history.");
    for (const auto& point : mesh.nodes)
        for (double value : point) if (!std::isfinite(value))
            throw std::invalid_argument("Non-finite mesh coordinate.");
    for (const auto& e : mesh.elements) for (auto n : e)
        if (n >= mesh.nodes.size()) throw std::invalid_argument("Invalid connectivity.");
    double previousTime = 0;
    bool first = true;
    for (const auto& record : result.history) {
        const auto& s = record.state;
        if (s.displacement.size()!=2*mesh.nodes.size() ||
            s.velocity.size()!=s.displacement.size() || s.acceleration.size()!=s.displacement.size() ||
            !std::isfinite(s.time) || (!first && s.time<=previousTime) ||
            !std::isfinite(record.kineticEnergy) || record.kineticEnergy<0 ||
            !std::isfinite(record.strainEnergy) || record.strainEnergy<0 ||
            !std::isfinite(record.kineticEnergy+record.strainEnergy))
            throw std::invalid_argument("Invalid dynamic record.");
        for (const auto* values : {&s.displacement,&s.velocity,&s.acceleration})
            for (double value : *values) if (!std::isfinite(value))
                throw std::invalid_argument("Non-finite dynamic state.");
        previousTime=s.time; first=false;
    }
    if (!std::filesystem::create_directories(directory))
        throw std::runtime_error("Output directory already exists.");
    auto nodes=open(directory/"nodes.csv");
    nodes<<"node,x,y\n";
    for(std::size_t n=0;n<mesh.nodes.size();++n)
        nodes<<n<<','<<mesh.nodes[n][0]<<','<<mesh.nodes[n][1]<<'\n';
    auto elements=open(directory/"elements.csv");
    elements<<"element,n0,n1,n2,n3\n";
    for(std::size_t i=0;i<mesh.elements.size();++i) {
        elements<<i;for(auto n:mesh.elements[i])elements<<','<<n;elements<<'\n';
    }
    auto history=open(directory/"history.csv");
    history<<"frame,time,node,ux,uy,vx,vy,ax,ay,kinetic_energy,strain_energy,total_energy\n";
    auto frames=open(directory/"frames.csv");
    frames<<"frame,time,node,ux,uy\n";
    for(std::size_t i=0;i<result.history.size();++i) {
        const auto& r=result.history[i];const auto& s=r.state;const auto n=2*observationNode;
        history<<i<<','<<s.time<<','<<observationNode<<','<<s.displacement[n]<<','<<s.displacement[n+1]
            <<','<<s.velocity[n]<<','<<s.velocity[n+1]<<','<<s.acceleration[n]<<','<<s.acceleration[n+1]
            <<','<<r.kineticEnergy<<','<<r.strainEnergy<<','<<r.kineticEnergy+r.strainEnergy<<'\n';
        for(std::size_t node=0;node<mesh.nodes.size();++node)
            frames<<i<<','<<s.time<<','<<node<<','<<s.displacement[2*node]<<','<<s.displacement[2*node+1]<<'\n';
    }
}
}
