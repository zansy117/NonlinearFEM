#include "nonlinearfem/Assembly.h"
#include "nonlinearfem/ExplicitDynamics.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>
using namespace NonlinearFEM;
namespace {
void require(bool condition,const char* message){if(!condition)throw std::runtime_error(message);}
struct Run {std::vector<DynamicVector> displacement;double energyError=0;};
Run simulate(int refinement){
    // Reference rectangle L=4, H=1; 4x2 Q4 elements, 15 shared nodes.
    Q4Mesh mesh{};
    for(int j=0;j<=2;++j)for(int i=0;i<=4;++i)mesh.nodes.push_back({double(i),.5*j});
    for(std::size_t j=0;j<2;++j)for(std::size_t i=0;i<4;++i){const auto a=5*j+i;mesh.elements.push_back({a,a+1,a+6,a+5});}
    MaterialParameters material{};material.youngModulus=1000;material.poissonRatio=.3;
    const double thickness=.2,density=1,amplitude=1e-4;
    const auto mass=assembleGlobalLumpedMass(mesh,density,thickness);
    const auto N=2*mesh.nodes.size();require(mass.size()==N,"Mass dimensions");
    double sumX=0,sumY=0;
    for(std::size_t n=0;n<mesh.nodes.size();++n){
        require(std::isfinite(mass[2*n])&&mass[2*n]>0&&mass[2*n]==mass[2*n+1],"Invalid assembled mass");
        sumX+=mass[2*n];sumY+=mass[2*n+1];
    }
    require(std::abs(sumX-.8)<1e-13&&std::abs(sumY-.8)<1e-13,"Reference mass conservation");
    std::vector<bool> fixed(N,false);
    ExplicitState state{};state.displacement.assign(N,0);state.velocity.assign(N,0);
    for(std::size_t n=0;n<mesh.nodes.size();++n){
        const double x=mesh.nodes[n][0],y=mesh.nodes[n][1]-.5;
        if(n%5==0){fixed[2*n]=fixed[2*n+1]=true;}
        // Smooth bending-like perturbation: w(0)=w'(0)=0.
        state.displacement[2*n]=-y*2*amplitude*x/16;
        state.displacement[2*n+1]=amplitude*x*x/16;
    }
    const InternalForceFunction internal=[&](const DynamicVector& u){return assembleInternalForce(mesh,u,material,thickness);};
    const ExternalForceFunction external=[N](double){return DynamicVector(N,0);};
    const auto energy=[&](const ExplicitState& s){
        double kinetic=0;for(std::size_t i=0;i<N;++i)kinetic+=.5*mass[i]*s.velocity[i]*s.velocity[i];
        return kinetic+assembleGlobalEnergy(mesh,s.displacement,material,thickness);
    };
    InitializeExplicitAcceleration(state,mass,fixed,internal,external);
    const double initial=energy(state);require(std::isfinite(initial)&&initial>0,"Initial energy");
    Run run;bool crossedZero=false;double maximum=0;
    const double dt=2e-4/refinement;
    const int steps=4000*refinement;
    for(int step=0;step<=steps;++step){
        require(std::abs(state.time-step*dt)<1e-10,"Time tracking");
        for(std::size_t i=0;i<N;++i){
            require(std::isfinite(state.displacement[i])&&std::isfinite(state.velocity[i])&&std::isfinite(state.acceleration[i]),"Non-finite dynamic state");
            if(fixed[i])require(state.displacement[i]==0&&state.velocity[i]==0&&state.acceleration[i]==0,"Fixed DOF drift");
            maximum=std::max(maximum,std::abs(state.displacement[i]));
        }
        const double current=energy(state);require(std::isfinite(current)&&current>=0,"Invalid total energy");
        run.energyError=std::max(run.energyError,std::abs(current-initial)/initial);
        crossedZero=crossedZero||state.displacement[19]<0; // Right-edge midpoint, node 9.
        if(step%(10*refinement)==0)run.displacement.push_back(state.displacement);
        if(step<steps)advanceCentralDifference(state,mass,fixed,dt,internal,external);
    }
    require(crossedZero,"No free-vibration sign reversal");
    require(maximum<5*amplitude,"Unbounded vibration response");
    require(run.displacement.size()==401,"Common sample count");
    require(run.energyError<1e-4,"Excessive total-energy drift");
    std::cout<<"dt="<<dt<<", maximum relative energy error="<<run.energyError<<'\n';
    return run;
}
double difference(const Run& a,const Run& b){
    require(a.displacement.size()==b.displacement.size(),"History dimensions");double error=0;
    for(std::size_t t=0;t<a.displacement.size();++t)for(std::size_t i=0;i<a.displacement[t].size();++i)
        error=std::max(error,std::abs(a.displacement[t][i]-b.displacement[t][i]));
    return error;
}
}
int main(){try{
    const auto a=simulate(1),b=simulate(2),c=simulate(4);
    const double ab=difference(a,b),bc=difference(b,c);
    std::cout<<"Common-time displacement differences="<<ab<<", "<<bc<<"; ratio="<<ab/bc<<'\n';
    std::cout<<"Energy refinement ratios="<<a.energyError/b.energyError<<", "<<b.energyError/c.energyError<<'\n';
    require(ab>1e-14&&bc>1e-14&&ab/bc>3.5&&ab/bc<4.5,"Displacement second-order refinement");
    require(bc/1e-4<1e-4,"Fine-step displacement change too large");
    require(a.energyError/b.energyError>3.3&&a.energyError/b.energyError<4.7&&b.energyError/c.energyError>3.3&&b.energyError/c.energyError<4.7,"Energy second-order refinement");
    std::cout<<"PASS: multi-element explicit free vibration, constraints, mass, bounded response, displacement and energy refinement.\n";
}catch(const std::exception& error){std::cerr<<error.what()<<'\n';return 1;}return 0;}