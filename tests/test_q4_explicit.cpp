#include "nonlinearfem/Assembly.h"
#include "nonlinearfem/ExplicitDynamics.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace NonlinearFEM;
namespace {
void require(bool condition,const char* message) {
    if(!condition) throw std::runtime_error(message);
}
using History=std::vector<DynamicVector>;
History simulate(int refinement) {
    const Q4Mesh mesh={{{0,0},{1,0},{1,1},{0,1}},{{0,1,2,3}}};
    MaterialParameters material{};
    material.youngModulus=1000;material.poissonRatio=.3;
    const double density=1,thickness=.2;
    const auto mass=assembleGlobalLumpedMass(mesh,density,thickness);
    require(mass.size()==8,"Mass dimensions");
    for(double m:mass) require(std::isfinite(m)&&std::abs(m-.05)<1e-14,"Unexpected unit-square mass");
    const std::vector<bool> fixed={true,true,false,false,false,false,true,true};
    ExplicitState state{};
    state.displacement.assign(8,0);state.velocity.assign(8,0);
    state.displacement[2]=state.displacement[4]=1e-4;
    const InternalForceFunction internal=[&](const DynamicVector& u) {
        return assembleInternalForce(mesh,u,material,thickness);
    };
    const ExternalForceFunction external=[](double){return DynamicVector(8,0);};
    InitializeExplicitAcceleration(state,mass,fixed,internal,external);
    require(state.acceleration[2]<0&&state.acceleration[4]<0,"Initial acceleration is not restoring");
    // Independent homogeneous initial-state force calculation: F=diag(1+eps,1), P=F*S.
    const double stretch=1.0001;
    const double exx=(stretch*stretch-1)/2;
    const double sxx=material.youngModulus*exx/(1-material.poissonRatio*material.poissonRatio);
    const double expectedAx=-stretch*sxx*thickness/(2*.05);
    require(std::abs(state.acceleration[2]-expectedAx)<1e-10,"Initial acceleration disagrees with analytic PK force");
    History history;
    bool crossedZero=false;
    const double dt=1e-4/refinement;
    const int steps=2000*refinement;
    double maxDisplacement=0;
    for(int step=0;step<=steps;++step) {
        require(std::abs(state.time-step*dt)<1e-11,"Incorrect physical time");
        for(std::size_t i=0;i<8;++i) {
            require(std::isfinite(state.displacement[i])&&std::isfinite(state.velocity[i])&&std::isfinite(state.acceleration[i]),"Non-finite dynamic state");
            if(fixed[i]) require(state.displacement[i]==0&&state.velocity[i]==0&&state.acceleration[i]==0,"Fixed DOF moved");
            maxDisplacement=std::max(maxDisplacement,std::abs(state.displacement[i]));
        }
        // Symmetry of this geometry and initial condition must persist.
        require(std::abs(state.displacement[2]-state.displacement[4])<1e-12,"Axial symmetry lost");
        require(std::abs(state.displacement[3]+state.displacement[5])<1e-12,"Transverse symmetry lost");
        crossedZero=crossedZero||state.displacement[2]<0;
        if(step%(10*refinement)==0) history.push_back(state.displacement);
        if(step<steps) advanceCentralDifference(state,mass,fixed,dt,internal,external);
    }
    require(crossedZero,"No free-vibration sign reversal");
    require(maxDisplacement<2e-4,"Unbounded small-amplitude response");
    require(history.size()==201,"Common sample times mismatch");
    return history;
}
double difference(const History& a,const History& b) {
    require(a.size()==b.size(),"History dimensions");
    double result=0;
    for(std::size_t t=0;t<a.size();++t)
        for(std::size_t i=0;i<a[t].size();++i)
            result=std::max(result,std::abs(a[t][i]-b[t][i]));
    return result;
}
}
int main() {
    try {
        const auto coarse=simulate(1),medium=simulate(2),fine=simulate(4);
        const double d1=difference(coarse,medium),d2=difference(medium,fine);
        std::cout<<"Common-time displacement differences: "<<d1<<", "<<d2<<"; ratio="<<d1/d2<<'\n';
        require(d1>1e-14&&d2>1e-14,"Refinement differences too small for a meaningful order check");
        require(d1/d2>3.5&&d1/d2<4.5,"Expected second-order temporal convergence");
        require(d2/1e-4<1e-4,"Fine-grid time-step difference too large");
        std::cout<<"PASS: Q4 mass/internal-force connection, analytic initial acceleration, constraints, symmetry, bounded oscillation and time-step refinement.\n";
    } catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 1;}
    return 0;
}