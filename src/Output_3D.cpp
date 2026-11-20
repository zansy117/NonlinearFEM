#include "nonlinearfem/ModelIO_3D.h"
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
void writeDynamicResults_3D(const H8Mesh& mesh,
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
        if (s.displacement.size()!=3*mesh.nodes.size() ||
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
    nodes<<"node,x,y,z\n";
    for(std::size_t n=0;n<mesh.nodes.size();++n)
        nodes<<n<<','<<mesh.nodes[n][0]<<','<<mesh.nodes[n][1]<<','<<mesh.nodes[n][2]<<'\n';
    auto elements=open(directory/"elements.csv");
    elements<<"element,n0,n1,n2,n3,n4,n5,n6,n7\n";
    for(std::size_t i=0;i<mesh.elements.size();++i) {
        elements<<i;for(auto n:mesh.elements[i])elements<<','<<n;elements<<'\n';
    }
    auto history=open(directory/"history.csv");
    history<<"frame,time,node,ux,uy,uz,vx,vy,vz,ax,ay,az,kinetic_energy,strain_energy,total_energy\n";
    auto frames=open(directory/"frames.csv");
    frames<<"frame,time,node,ux,uy,uz\n";
    for(std::size_t i=0;i<result.history.size();++i) {
        const auto& r=result.history[i];const auto& s=r.state;const auto n=3*observationNode;
        history<<i<<','<<s.time<<','<<observationNode<<','<<s.displacement[n]<<','<<s.displacement[n+1]<<','<<s.displacement[n+2]
            <<','<<s.velocity[n]<<','<<s.velocity[n+1]<<','<<s.velocity[n+2]<<','<<s.acceleration[n]<<','<<s.acceleration[n+1]<<','<<s.acceleration[n+2]
            <<','<<r.kineticEnergy<<','<<r.strainEnergy<<','<<r.kineticEnergy+r.strainEnergy<<'\n';
        for(std::size_t node=0;node<mesh.nodes.size();++node)
            frames<<i<<','<<s.time<<','<<node<<','<<s.displacement[3*node]<<','<<s.displacement[3*node+1]<<','<<s.displacement[3*node+2]<<'\n';
    }
}
}

namespace NonlinearFEM {
void writeResults_3D(const AnalysisModel_3D& model,const LoadSteppingResult_3D& result,
    const std::filesystem::path& directory) {
    if(!result.converged || std::abs(result.finalLoadFactor-1)>1e-10)
        throw std::runtime_error("Only a fully converged H8 static state can be exported.");
    DynamicAnalysisResult snapshots;
    DynamicRecord record;record.state.displacement=result.displacement;
    record.state.velocity.assign(result.displacement.size(),0);
    record.state.acceleration.assign(result.displacement.size(),0);
    record.strainEnergy=assembleGlobalEnergy_3D(model.mesh,result.displacement,model.material);
    snapshots.history.push_back(record);
    writeDynamicResults_3D(model.mesh,snapshots,0,directory);
    const auto internal=assembleInternalForce_3D(model.mesh,result.displacement,model.material);
    const auto external=assembleExternalForce_3D(model.mesh.nodes.size(),model.loads);
    std::vector<bool> fixed(internal.size(),false);
    for(const auto& c:model.constraints)fixed[3*c.node+c.direction]=true;
    auto file=open(directory/"reactions.csv");file<<"node,rx,ry,rz\n";
    for(std::size_t n=0;n<model.mesh.nodes.size();++n){file<<n;
        for(std::size_t d=0;d<3;++d){auto i=3*n+d;file<<','<<(fixed[i]?internal[i]-external[i]:0);}
        file<<'\n';}
    auto steps=open(directory/"load_history.csv");steps<<"step,load_factor,increment,iterations,retries,residual\n";
    for(std::size_t i=0;i<result.history.size();++i){const auto& h=result.history[i];
        steps<<i+1<<','<<h.loadFactor<<','<<h.loadIncrement<<','<<h.iterations<<','<<h.retries<<','<<h.residualNorm<<'\n';}
    auto info=open(directory/"analysis.txt");info<<"analysis static\ndimension 3\nload_factor "<<result.finalLoadFactor<<"\n";
}
}
